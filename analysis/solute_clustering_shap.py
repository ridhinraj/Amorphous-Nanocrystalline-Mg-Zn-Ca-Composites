"""
solute_clustering_shap.py

Identifies which local Zn/Ca solute configurations are statistically
associated with high-stress (shear-band) regions, using per-atom
features and an explainable gradient-boosting classifier.

------------------------------------------------------------------
DESIGN NOTE -- what this uses as a label, and why

The original plan called for labeling atoms by proximity to a
PTM-identified crystalline/amorphous interphase boundary. That
requires a working PTM phase-identification step, which hasn't been
built in this project yet (OVITO's Polyhedral Template Matching was
planned but never implemented against real data). Rather than block
this script on a piece of tooling that doesn't exist, it reuses what
detect_shear_bands.py already computes: per-atom von Mises stress
from the tension dump. Atoms are labeled "near-band" if they sit
within a distance threshold of a high-von-Mises ("hot") atom -- a
legitimate proxy for "near an active zone of localized deformation",
though not identical to "near the amorphous/crystalline interphase
specifically". This should be revisited once real PTM identification
exists; until then, treat "near-band" as "near a high-stress
concentration", which is what it actually measures.

FEATURES computed per Zn/Ca atom:
  - n_Zn_neighbors, n_Ca_neighbors, n_Mg_neighbors within a cutoff
    (3.25 Angstrom, chosen to match the first-neighbor-shell distance
    established from this project's RDF checks)
  - own_type (Zn or Ca)
  - distance to the single nearest high-von-Mises atom (used to
    derive the label, then EXCLUDED from the feature set fed to the
    classifier -- including it would let the model trivially recover
    the label instead of learning from clustering chemistry)

TARGET: binary, "near-band" (1) vs "far-from-band" (0).
------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import shap

from detect_shear_bands import parse_dump, von_mises_stress

NEIGHBOR_CUTOFF = 3.25   # Angstrom, first-shell distance from this project's RDF checks
BAND_PROXIMITY_CUTOFF = 6.0  # Angstrom -- "near" a hot atom
HOT_PERCENTILE = 90.0


def compute_solute_features(atoms_df, neighbor_cutoff=NEIGHBOR_CUTOFF,
                             band_proximity_cutoff=BAND_PROXIMITY_CUTOFF,
                             hot_percentile=HOT_PERCENTILE):
    """
    Given one frame's atoms_df (columns: id, type, x, y, z, sxx..syz,
    matching detect_shear_bands.parse_dump's output), returns a
    DataFrame of per-Zn/Ca-atom features and the near-band label.

    Type mapping (project-wide convention): 1=Mg, 2=Zn, 3=Ca.
    """
    positions = atoms_df[["x", "y", "z"]].values
    types = atoms_df["type"].values

    vm = von_mises_stress(
        atoms_df["sxx"], atoms_df["syy"], atoms_df["szz"],
        atoms_df["sxy"], atoms_df["sxz"], atoms_df["syz"],
    ).values
    hot_threshold = np.percentile(vm, hot_percentile)
    hot_positions = positions[vm >= hot_threshold]

    tree_all = cKDTree(positions)
    tree_hot = cKDTree(hot_positions) if len(hot_positions) > 0 else None

    solute_mask = (types == 2) | (types == 3)  # Zn or Ca
    solute_idx = np.where(solute_mask)[0]

    records = []
    for i in solute_idx:
        pos = positions[i]
        neighbor_idx = tree_all.query_ball_point(pos, r=neighbor_cutoff)
        neighbor_idx = [j for j in neighbor_idx if j != i]
        neighbor_types = types[neighbor_idx]

        if tree_hot is not None and len(hot_positions) > 0:
            dist_to_hot, _ = tree_hot.query(pos)
        else:
            dist_to_hot = np.inf

        records.append({
            "atom_id": int(atoms_df.iloc[i]["id"]),
            "own_type_Ca": 1 if types[i] == 3 else 0,  # 1=Ca, 0=Zn
            "n_Mg_neighbors": int(np.sum(neighbor_types == 1)),
            "n_Zn_neighbors": int(np.sum(neighbor_types == 2)),
            "n_Ca_neighbors": int(np.sum(neighbor_types == 3)),
            "dist_to_hot_atom": float(dist_to_hot),
            "near_band": int(dist_to_hot <= band_proximity_cutoff),
        })

    return pd.DataFrame(records)


def train_and_explain(features_df, feature_columns=None, test_size=0.25, random_state=0):
    """
    Trains a gradient-boosting classifier on near_band, holds out a
    test split for an honest AUC check, then computes SHAP values on
    the test set to identify which features actually drive the
    near-band classification.
    """
    if feature_columns is None:
        feature_columns = ["own_type_Ca", "n_Mg_neighbors", "n_Zn_neighbors", "n_Ca_neighbors"]

    X = features_df[feature_columns]
    y = features_df["near_band"]

    if y.nunique() < 2:
        raise ValueError(
            "near_band has only one class in this data -- can't train a "
            "classifier. Check band_proximity_cutoff/hot_percentile, or "
            "this frame genuinely has no meaningful high-stress region yet."
        )

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    clf = GradientBoostingClassifier(random_state=random_state)
    clf.fit(X_train, y_train)

    test_auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])

    explainer = shap.TreeExplainer(clf)
    shap_values = explainer.shap_values(X_test)

    mean_abs_shap = pd.Series(
        np.abs(shap_values).mean(axis=0), index=feature_columns
    ).sort_values(ascending=False)

    return clf, test_auc, mean_abs_shap, (X_test, shap_values)


def compute_solute_features_continuous(atoms_df, neighbor_cutoff=NEIGHBOR_CUTOFF):
    """
    Regression variant, added after the classification version's
    "near_band" label turned out to be near-degenerate on real data
    (97.4% positive rate on composition 5 -- consistent with
    detect_shear_bands.py's finding that this system's stress field
    is diffuse rather than spatially localized into discrete bands).
    A binary "near/far from a band" label needs genuine spatial
    localization to be meaningful; this system doesn't have it.

    Instead of a binary label, this predicts each solute atom's OWN
    von Mises stress directly (log-transformed, since stress-volume
    values span orders of magnitude) from its local Zn/Ca/Mg
    coordination. This doesn't assume band-like structure at all --
    it asks whether local solute clustering correlates with local
    stress magnitude, atom by atom, which is a valid question
    regardless of whether discrete bands exist. It also connects
    directly to the bulk-scale Ca/yield-stress trend already found,
    at the atomic scale instead.
    """
    positions = atoms_df[["x", "y", "z"]].values
    types = atoms_df["type"].values
    vm = von_mises_stress(
        atoms_df["sxx"], atoms_df["syy"], atoms_df["szz"],
        atoms_df["sxy"], atoms_df["sxz"], atoms_df["syz"],
    ).values

    tree_all = cKDTree(positions)
    solute_mask = (types == 2) | (types == 3)
    solute_idx = np.where(solute_mask)[0]

    records = []
    for i in solute_idx:
        pos = positions[i]
        neighbor_idx = tree_all.query_ball_point(pos, r=neighbor_cutoff)
        neighbor_idx = [j for j in neighbor_idx if j != i]
        neighbor_types = types[neighbor_idx]

        records.append({
            "atom_id": int(atoms_df.iloc[i]["id"]),
            "own_type_Ca": 1 if types[i] == 3 else 0,
            "n_Mg_neighbors": int(np.sum(neighbor_types == 1)),
            "n_Zn_neighbors": int(np.sum(neighbor_types == 2)),
            "n_Ca_neighbors": int(np.sum(neighbor_types == 3)),
            "log_von_mises": float(np.log(max(vm[i], 1.0))),
        })

    return pd.DataFrame(records)


def train_and_explain_regression(features_df, feature_columns=None, test_size=0.25, random_state=0):
    """
    Same SHAP-explainability approach as train_and_explain, but for
    the continuous log_von_mises target instead of a binary label.
    Reports held-out R^2 instead of AUC.
    """
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.metrics import r2_score

    if feature_columns is None:
        feature_columns = ["own_type_Ca", "n_Mg_neighbors", "n_Zn_neighbors", "n_Ca_neighbors"]

    X = features_df[feature_columns]
    y = features_df["log_von_mises"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )

    reg = GradientBoostingRegressor(random_state=random_state)
    reg.fit(X_train, y_train)
    test_r2 = r2_score(y_test, reg.predict(X_test))

    explainer = shap.TreeExplainer(reg)
    shap_values = explainer.shap_values(X_test)
    mean_abs_shap = pd.Series(
        np.abs(shap_values).mean(axis=0), index=feature_columns
    ).sort_values(ascending=False)

    return reg, test_r2, mean_abs_shap, (X_test, shap_values)


# ============================================================
# Self-test against synthetic data with a KNOWN, engineered
# correlation: Ca atoms are deliberately placed preferentially
# near the synthetic "hot" region, so a working pipeline should
# recover own_type_Ca and/or n_Ca_neighbors as important features.
# ============================================================
def _write_synthetic_dump_with_solute_pattern(path, n_atoms=6000, box_len=None, seed=0):
    if box_len is None:
        real_density = 389344 / (207.0 ** 3)
        box_len = (n_atoms / real_density) ** (1.0 / 3.0)

    rng = np.random.default_rng(seed)
    positions = rng.uniform(0, box_len, size=(n_atoms, 3))

    # Hot region: a slab near x=0, same construction as detect_shear_bands' test.
    hot_region_mask = positions[:, 0] < 0.15 * box_len

    types = np.full(n_atoms, 1, dtype=int)  # default Mg
    n_zn = int(0.30 * n_atoms)
    n_ca = int(0.05 * n_atoms)

    # Ca preferentially placed IN the hot region (engineered correlation);
    # Zn placed uniformly at random, no correlation with the hot region.
    ca_in_hot = rng.choice(np.where(hot_region_mask)[0],
                            size=min(n_ca, int(hot_region_mask.sum())), replace=False)
    remaining_ca = n_ca - len(ca_in_hot)
    other_idx = np.where(~np.isin(np.arange(n_atoms), ca_in_hot))[0]
    ca_elsewhere = rng.choice(other_idx, size=max(remaining_ca, 0), replace=False)
    ca_idx = np.concatenate([ca_in_hot, ca_elsewhere])
    types[ca_idx] = 3

    remaining_idx = np.setdiff1d(np.arange(n_atoms), ca_idx)
    zn_idx = rng.choice(remaining_idx, size=n_zn, replace=False)
    types[zn_idx] = 2

    stress = rng.normal(0, 1.0, size=(n_atoms, 6))
    stress[hot_region_mask] = rng.normal(20.0, 1.0, size=(hot_region_mask.sum(), 6))

    timestep = 20000
    with open(path, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{timestep}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{n_atoms}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write(f"0.0 {box_len}\n0.0 {box_len}\n0.0 {box_len}\n")
        f.write("ITEM: ATOMS id type x y z sxx syy szz sxy sxz syz\n")
        for i in range(n_atoms):
            x, y, z = positions[i]
            sxx, syy, szz, sxy, sxz, syz = stress[i]
            f.write(f"{i+1} {types[i]} {x:.4f} {y:.4f} {z:.4f} "
                    f"{sxx:.4f} {syy:.4f} {szz:.4f} {sxy:.4f} {sxz:.4f} {syz:.4f}\n")


def _self_test():
    print("=" * 60)
    print("SELF-TEST against synthetic dump -- not real results")
    print("=" * 60)

    synth_path = "synthetic_solute_dump.dump"
    _write_synthetic_dump_with_solute_pattern(synth_path)
    print(f"Wrote synthetic dump: {synth_path}")

    frames = list(parse_dump(synth_path))
    timestep, box_bounds, atoms_df = frames[0]

    features_df = compute_solute_features(atoms_df)
    print(f"\nComputed features for {len(features_df)} Zn/Ca atoms")
    print(f"near_band positive rate: {features_df['near_band'].mean():.3f}")

    clf, test_auc, mean_abs_shap, _ = train_and_explain(features_df)

    print(f"\nHeld-out test AUC: {test_auc:.3f}")
    print("\nMean |SHAP value| per feature (higher = more important):")
    print(mean_abs_shap.to_string())

    print()
    top_feature = mean_abs_shap.index[0]
    if top_feature in ("own_type_Ca", "n_Ca_neighbors"):
        print(f"PASS: top feature is '{top_feature}', matching the engineered "
              f"Ca-near-hot-region correlation.")
    else:
        print(f"WARNING: top feature is '{top_feature}', not a Ca-related feature -- "
              f"expected 'own_type_Ca' or 'n_Ca_neighbors' given how this test was built. "
              f"Check band_proximity_cutoff/hot_percentile before trusting this on real data.")

    assert test_auc > 0.6, f"Expected clearly-better-than-random AUC, got {test_auc:.3f}"
    print("\nAssertion passed (AUC > 0.6).")


# ============================================================
# NEXT STEP: once a real tension_trajectory.dump exists, run:
#
#   from detect_shear_bands import parse_dump
#   from solute_clustering_shap import compute_solute_features, train_and_explain
#
#   frames = list(parse_dump("tension_trajectory.dump"))
#   timestep, box_bounds, atoms_df = frames[-1]   # last (most strained) frame
#   features_df = compute_solute_features(atoms_df)
#   clf, auc, shap_importance, _ = train_and_explain(features_df)
#
# Same caveat as detect_shear_bands.py: NEIGHBOR_CUTOFF, BAND_PROXIMITY_
# CUTOFF, and HOT_PERCENTILE are reasonable starting points given this
# project's real atomic spacing, not yet tuned against real data.
# ============================================================

# ============================================================
# Self-test for the regression variant: engineered CONTINUOUS
# correlation between local Ca clustering and stress magnitude
# (unlike the classification self-test's binary hot-region split),
# since that's the actual shape of the question now being asked.
# ============================================================
def _write_synthetic_dump_with_continuous_correlation(path, n_atoms=6000, box_len=None, seed=0):
    if box_len is None:
        real_density = 389344 / (207.0 ** 3)
        box_len = (n_atoms / real_density) ** (1.0 / 3.0)

    rng = np.random.default_rng(seed)
    positions = rng.uniform(0, box_len, size=(n_atoms, 3))

    types = np.full(n_atoms, 1, dtype=int)
    n_zn = int(0.30 * n_atoms)
    n_ca = int(0.05 * n_atoms)
    idx = rng.permutation(n_atoms)
    types[idx[:n_ca]] = 3
    types[idx[n_ca:n_ca + n_zn]] = 2

    # Engineer: atoms with more Ca neighbors get systematically LOWER
    # stress (a stand-in for "Ca clustering locally relieves stress",
    # a plausible mechanism consistent with the real Ca/yield trend).
    tree = cKDTree(positions)
    stress = rng.normal(0, 1.0, size=(n_atoms, 6))
    for i in range(n_atoms):
        if types[i] in (2, 3):
            neighbor_idx = tree.query_ball_point(positions[i], r=NEIGHBOR_CUTOFF)
            n_ca_nearby = np.sum(types[neighbor_idx] == 3)
            reduction = 2.0 * n_ca_nearby  # more Ca neighbors -> lower stress
            stress[i] = rng.normal(10.0 - reduction, 1.0, size=6)

    timestep = 20000
    with open(path, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{timestep}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{n_atoms}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write(f"0.0 {box_len}\n0.0 {box_len}\n0.0 {box_len}\n")
        f.write("ITEM: ATOMS id type x y z sxx syy szz sxy sxz syz\n")
        for i in range(n_atoms):
            x, y, z = positions[i]
            sxx, syy, szz, sxy, sxz, syz = stress[i]
            f.write(f"{i+1} {types[i]} {x:.4f} {y:.4f} {z:.4f} "
                    f"{sxx:.4f} {syy:.4f} {szz:.4f} {sxy:.4f} {sxz:.4f} {syz:.4f}\n")


def _self_test_regression():
    print("=" * 60)
    print("SELF-TEST (regression variant) -- not real results")
    print("=" * 60)

    synth_path = "synthetic_solute_regression_dump.dump"
    _write_synthetic_dump_with_continuous_correlation(synth_path)
    print(f"Wrote synthetic dump: {synth_path}")

    frames = list(parse_dump(synth_path))
    timestep, box_bounds, atoms_df = frames[0]

    features_df = compute_solute_features_continuous(atoms_df)
    print(f"\nComputed features for {len(features_df)} Zn/Ca atoms")
    print(f"log_von_mises range: {features_df['log_von_mises'].min():.2f} to "
          f"{features_df['log_von_mises'].max():.2f}")

    reg, test_r2, shap_importance, _ = train_and_explain_regression(features_df)

    print(f"\nHeld-out test R^2: {test_r2:.3f}")
    print("\nMean |SHAP value| per feature:")
    print(shap_importance.to_string())

    top_feature = shap_importance.index[0]
    print()
    if top_feature == "n_Ca_neighbors":
        print(f"PASS: top feature is '{top_feature}', matching the engineered "
              f"Ca-neighbor-count vs. stress correlation.")
    else:
        print(f"WARNING: top feature is '{top_feature}', not 'n_Ca_neighbors' -- "
              f"expected given how this test was built.")

    assert test_r2 > 0.3, f"Expected clearly-better-than-baseline R^2, got {test_r2:.3f}"
    print("\nAssertion passed (R^2 > 0.3).")


if __name__ == "__main__":
    _self_test()
    print()
    _self_test_regression()
