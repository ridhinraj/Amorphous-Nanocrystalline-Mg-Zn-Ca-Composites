"""
gp_active_learning.py

Gaussian process surrogate + active learning for the AN composite
mechanical-property dataset. Built and tested against SYNTHETIC data
(see `_synthetic_demo` at the bottom) since round-1 MD data doesn't
exist yet -- point `load_training_table` at the real CSV the moment
it does, and everything else here needs no changes.

------------------------------------------------------------------
TRAINING TABLE SCHEMA (locked now, before real data exists, so the
columns from round 1 and round 2 stay consistent):

    run_id            str    e.g. "run_1", "vf_2", "round2_003"
    xZn, xCa, xMg      float  composition, mole fraction
    volume_fraction    float  0.0-1.0 (0.0 for the 10 composition-
                               sweep runs, which are ~0.7% VF single
                               seeds -- keep the real small number,
                               don't round to 0, so the GP still
                               sees the true input)
    seed_radius        float  Angstrom (composition-sweep runs) or
                               mean seed radius (VF-sweep runs)
    temperature_K      float
    strain_rate        float  1/s
    loading_mode       str    "tension" or "shear"
    yield_stress_GPa   float  TARGET
    band_count         float  TARGET (post-yield shear band count)
    band_spacing_nm    float  TARGET (mean spacing, NaN if band_count<2)
    arrest_strain      float  TARGET (strain at first band arrest)
    round               int   1 or 2 (which active-learning round)
    source             str    "fresh_quench" or "reused_glass"
    job_id             str    Slurm job ID, for provenance
------------------------------------------------------------------
"""

import warnings
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from sklearn.preprocessing import StandardScaler

FEATURE_COLUMNS = [
    "xZn", "xCa", "volume_fraction", "seed_radius",
    "temperature_K", "strain_rate", "loading_mode_shear",
]
TARGET_COLUMNS = ["yield_stress_GPa", "band_count", "band_spacing_nm", "arrest_strain"]


def _encode_loading_mode(df):
    """Shared by every code path that produces a training/candidate
    table (real CSV load, synthetic demo, candidate pool) so the
    encoding can't silently diverge between them."""
    df = df.copy()
    df["loading_mode_shear"] = (df["loading_mode"] == "shear").astype(float)
    return df


def load_training_table(csv_path):
    """Load the training table and encode loading_mode as 0/1."""
    df = pd.read_csv(csv_path)
    return _encode_loading_mode(df)


def build_gp(length_scale_init=1.0, noise_init=1e-3):
    """
    Matern kernel (nu=2.5) is a reasonable default for physical
    response surfaces -- smoother than nu=1.5, doesn't assume the
    infinite differentiability of an RBF kernel (which tends to
    over-smooth sharp transitions like a percolation threshold).
    WhiteKernel accounts for MD's inherent stochastic noise
    (thermostat/seed dependence) rather than treating every point
    as an exact observation.
    """
    kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(
        length_scale=length_scale_init, length_scale_bounds=(1e-2, 1e2), nu=2.5
    ) + WhiteKernel(noise_level=noise_init, noise_level_bounds=(1e-6, 1e1))
    return GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=5)


def fit_target(df, target_col, feature_scaler=None):
    """
    Fit one GP for one target column. Features are standardized
    (mean 0, std 1) since the raw scales are wildly different
    (strain_rate ~1e9, volume_fraction ~0-0.15) and an unscaled
    Matern kernel would be dominated by whichever feature happens
    to have the largest raw numbers.
    """
    X = df[FEATURE_COLUMNS].values
    y = df[target_col].values
    mask = ~np.isnan(y)  # band_spacing_nm is NaN when band_count < 2
    X, y = X[mask], y[mask]

    if feature_scaler is None:
        feature_scaler = StandardScaler().fit(X)
    Xs = feature_scaler.transform(X)

    gp = build_gp()
    with warnings.catch_warnings():
        # Benign on very clean data (including this synthetic demo,
        # by construction): the optimizer pushes the noise term to
        # its lower bound when the true noise really is near zero.
        # Real MD data has genuine thermostat/seed stochasticity, so
        # this is less likely to fire there -- but if it does, it's
        # worth checking noise_level_bounds isn't clipping a real
        # non-trivial noise estimate.
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        gp.fit(Xs, y)
    return gp, feature_scaler


def predict(gp, feature_scaler, candidates_df):
    X = candidates_df[FEATURE_COLUMNS].values
    Xs = feature_scaler.transform(X)
    mean, std = gp.predict(Xs, return_std=True)
    return mean, std


def propose_next_batch(gp, feature_scaler, candidate_pool_df, n_select=3):
    """
    Uncertainty-sampling acquisition: rank candidates by predictive
    std, take the top n_select. This is deliberately simple --
    greedy variance maximization rather than a joint-diversity
    criterion (e.g. via the predictive covariance matrix). That's
    a reasonable, defensible starting point and worth noting as
    such in the Methods section; upgrading to a batch-diversity
    criterion is a natural reviewer suggestion to have an answer
    ready for, not a requirement for this paper.
    """
    mean, std = predict(gp, feature_scaler, candidate_pool_df)
    order = np.argsort(-std)  # descending uncertainty
    selected = candidate_pool_df.iloc[order[:n_select]].copy()
    selected["predicted_mean"] = mean[order[:n_select]]
    selected["predicted_std"] = std[order[:n_select]]
    return selected


def loo_calibration_check(df, target_col):
    """
    Leave-one-out calibration check: for each training point, fit
    on all others, predict at the held-out point, standardize the
    residual by the predicted std. If the GP's uncertainty is
    honest, these standardized residuals should look roughly like
    a standard normal (~68% within +-1, ~95% within +-1.96).

    With only 13 round-1 points this will be a noisy diagnostic --
    report it anyway, honestly, rather than skip it. A GP that
    looks poorly calibrated on 13 points is a real limitation to
    state, not a reason to hide the check.
    """
    X_full = df[FEATURE_COLUMNS].values
    y_full = df[target_col].values
    mask = ~np.isnan(y_full)
    X_full, y_full = X_full[mask], y_full[mask]
    n = len(y_full)

    standardized_residuals = np.zeros(n)
    for i in range(n):
        train_idx = np.array([j for j in range(n) if j != i])
        scaler = StandardScaler().fit(X_full[train_idx])
        Xs_train = scaler.transform(X_full[train_idx])
        Xs_test = scaler.transform(X_full[i : i + 1])

        gp = build_gp()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=ConvergenceWarning)
            gp.fit(Xs_train, y_full[train_idx])
        mean, std = gp.predict(Xs_test, return_std=True)
        std = max(std[0], 1e-8)
        standardized_residuals[i] = (y_full[i] - mean[0]) / std

    frac_within_1sigma = np.mean(np.abs(standardized_residuals) <= 1.0)
    frac_within_2sigma = np.mean(np.abs(standardized_residuals) <= 1.96)
    print(f"LOO calibration for {target_col} (n={n}):")
    print(f"  standardized residual mean={standardized_residuals.mean():.3f}, "
          f"std={standardized_residuals.std():.3f} (ideal: 0, 1)")
    print(f"  fraction within +-1 sigma: {frac_within_1sigma:.2f} (ideal: ~0.68)")
    print(f"  fraction within +-1.96 sigma: {frac_within_2sigma:.2f} (ideal: ~0.95)")
    return standardized_residuals


def estimate_active_learning_savings(round1_n_compositions, round2_n_reused,
                                      days_per_composition=14, days_per_reused=1):
    """
    The concrete compute-savings number for Table 4 / the paper's
    data-efficiency claim: what round 2 actually cost, vs. what it
    would have cost if each round-2 point had instead been a fresh
    composition (the naive alternative).
    """
    actual_cost = round2_n_reused * days_per_reused
    naive_cost = round2_n_reused * days_per_composition
    savings_pct = 100.0 * (naive_cost - actual_cost) / naive_cost
    print(f"Round 2: {round2_n_reused} reused-glass configurations")
    print(f"  Actual cost: {actual_cost} days")
    print(f"  Cost if each were a fresh composition instead: {naive_cost} days")
    print(f"  Savings: {savings_pct:.1f}%")
    return savings_pct


# ============================================================
# Self-test against synthetic data. Run this file directly to
# confirm the pipeline works end-to-end before pointing it at
# real data.
# ============================================================
def _make_synthetic_training_table(n=13, seed=0):
    """
    Fabricated data with a KNOWN functional form, so the GP's
    predictions can be checked against ground truth -- this is
    purely for testing the pipeline, not a stand-in for real
    results anywhere in the paper.
    """
    rng = np.random.default_rng(seed)
    rows = []
    compositions = [
        (0.35, 0.05), (0.30, 0.04), (0.40, 0.06), (0.25, 0.05), (0.35, 0.03),
        (0.35, 0.07), (0.30, 0.06), (0.40, 0.04), (0.32, 0.05), (0.38, 0.05),
    ]
    for i, (xZn, xCa) in enumerate(compositions):
        vf = 0.007  # ~single-seed VF from the composition sweep
        radius = rng.uniform(15, 35)
        true_yield = 0.5 + 2.0 * xZn - 1.0 * xCa + 0.1 * rng.standard_normal()
        rows.append(dict(
            run_id=f"run_{i+1}", xZn=xZn, xCa=xCa, xMg=1 - xZn - xCa,
            volume_fraction=vf, seed_radius=radius, temperature_K=300.0,
            strain_rate=1e9, loading_mode="tension",
            yield_stress_GPa=true_yield,
            band_count=1 + rng.integers(0, 2),
            band_spacing_nm=np.nan, arrest_strain=0.05 + 0.01 * rng.standard_normal(),
            round=1, source="fresh_quench", job_id=f"job_{1000+i}",
        ))
    # 3 VF-sweep points on the baseline composition (run_1's xZn/xCa)
    for i, vf in enumerate([0.05, 0.10, 0.15]):
        true_yield = 0.5 + 2.0 * 0.35 - 1.0 * 0.05 + 3.0 * vf + 0.1 * rng.standard_normal()
        rows.append(dict(
            run_id=f"vf_{i+1}", xZn=0.35, xCa=0.05, xMg=0.60,
            volume_fraction=vf, seed_radius=12.0, temperature_K=300.0,
            strain_rate=1e9, loading_mode="tension",
            yield_stress_GPa=true_yield,
            band_count=max(1, int(4 * (1 - vf / 0.15))),
            band_spacing_nm=8.0 + 20.0 * vf + rng.standard_normal(),
            arrest_strain=0.05 + 0.02 * vf + 0.01 * rng.standard_normal(),
            round=1, source="reused_glass", job_id=f"job_{2000+i}",
        ))
    return _encode_loading_mode(pd.DataFrame(rows))


def _make_candidate_pool(baseline_xZn=0.35, baseline_xCa=0.05):
    """Candidate round-2 points: new VF/T/rate on the baseline composition only."""
    rows = []
    for vf in [0.02, 0.08, 0.12, 0.20, 0.25]:
        for T in [300.0, 310.0]:
            for rate in [1e8, 1e9]:
                rows.append(dict(
                    xZn=baseline_xZn, xCa=baseline_xCa,
                    volume_fraction=vf, seed_radius=12.0,
                    temperature_K=T, strain_rate=rate, loading_mode="tension",
                ))
    return _encode_loading_mode(pd.DataFrame(rows))


def _synthetic_demo():
    print("=" * 60)
    print("SYNTHETIC SELF-TEST -- not real results, pipeline check only")
    print("=" * 60)

    df = _make_synthetic_training_table()
    print(f"\nRound 1 training table: {len(df)} rows")
    print(df[["run_id", "xZn", "xCa", "volume_fraction", "yield_stress_GPa"]].to_string(index=False))

    gp, scaler = fit_target(df, "yield_stress_GPa")
    print("\nFitted GP kernel:", gp.kernel_)

    print()
    loo_calibration_check(df, "yield_stress_GPa")

    candidates = _make_candidate_pool()
    print(f"\nCandidate pool for round 2: {len(candidates)} configurations")
    selected = propose_next_batch(gp, scaler, candidates, n_select=3)
    print("\nTop 3 proposed round-2 configurations (highest predictive uncertainty):")
    print(selected[["volume_fraction", "temperature_K", "strain_rate",
                     "predicted_mean", "predicted_std"]].to_string(index=False))

    print()
    estimate_active_learning_savings(round1_n_compositions=10, round2_n_reused=3)


if __name__ == "__main__":
    _synthetic_demo()
