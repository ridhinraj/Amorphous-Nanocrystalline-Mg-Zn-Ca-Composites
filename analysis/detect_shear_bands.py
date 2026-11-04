"""
detect_shear_bands.py

Automated shear-band detection from LAMMPS tension-test dump files,
replacing manual measurement.

------------------------------------------------------------------
WHAT THIS ACTUALLY USES, AND WHY

week4_tension.in's dump command is:
    dump 1 all custom 5000 tension_trajectory.dump id type x y z \
        c_peratom_stress[1] c_peratom_stress[2] c_peratom_stress[3] \
        c_peratom_stress[4] c_peratom_stress[5] c_peratom_stress[6]

That's per-atom STRESS (sxx, syy, szz, sxy, sxz, syz), not strain --
there's no reference-configuration displacement-gradient calculation
anywhere in the production pipeline (that would need something like
OVITO's Atomic Strain modifier comparing every frame back to frame 0,
which was never set up). So this script computes per-atom von Mises
STRESS directly from what's actually in the dump, rather than
assuming a strain field that doesn't exist in the data. High-von-
Mises-stress regions are the atomistic signature of localized shear
in these simulations, so this is a legitimate (and simpler, more
honest) substitute for strain-based detection given what's on disk.

UNIT CAVEAT: LAMMPS's compute stress/atom outputs stress*volume
(energy units), not raw stress -- dividing by each atom's volume to
get calibrated GPa would need per-atom Voronoi volumes, which aren't
computed anywhere in this pipeline. This script works in raw
stress-volume units throughout, which is fine for percentile-based
thresholding and clustering (relative comparison, not absolute
values). A rough GPa conversion using box-volume / atom-count as a
global average atomic volume is provided as an optional, clearly
approximate convenience -- don't treat it as calibrated per-atom
stress.

STRAIN-TIME MAPPING: week4_tension.in resets the timestep to 0,
equilibrates for 20000 steps, THEN freezes L0 and starts the
erate=0.001 (1/ps) deformation with dt=0.001 ps. So strain=0
corresponds to timestep 20000, not 0 -- frame_strain = (timestep -
20000) * erate * dt. This is parameterized below (T0_STEPS,
ERATE_PER_PS, TIMESTEP_PS) rather than hardcoded silently, since
getting this wrong would misalign every band-vs-strain plot.
------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# Must match week4_tension.in exactly -- see docstring above.
T0_STEPS = 20000        # timestep at which L0 was frozen (strain=0)
ERATE_PER_PS = 0.001     # fix deform erate, 1/ps
TIMESTEP_PS = 0.001      # LAMMPS timestep, ps


def frame_strain(timestep, erate_per_ps=ERATE_PER_PS):
    return (timestep - T0_STEPS) * erate_per_ps * TIMESTEP_PS


def von_mises_stress(sxx, syy, szz, sxy, sxz, syz):
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2)
    )


def parse_dump(path):
    """
    Yields (timestep, box_bounds, atoms_df) for each frame in a
    standard LAMMPS custom dump file matching week4_tension.in's
    column order: id type x y z sxx syy szz sxy sxz syz.
    """
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(f"Unexpected dump format, expected ITEM: TIMESTEP, got: {line!r}")
            timestep = int(f.readline().strip())

            f.readline()  # ITEM: NUMBER OF ATOMS
            n_atoms = int(f.readline().strip())

            f.readline()  # ITEM: BOX BOUNDS ...
            box_bounds = np.array([
                [float(v) for v in f.readline().split()] for _ in range(3)
            ])

            f.readline()  # ITEM: ATOMS ...

            rows = [f.readline().split() for _ in range(n_atoms)]
            arr = np.array(rows, dtype=float)
            atoms_df = pd.DataFrame(arr, columns=[
                "id", "type", "x", "y", "z",
                "sxx", "syy", "szz", "sxy", "sxz", "syz",
            ])
            atoms_df["id"] = atoms_df["id"].astype(int)
            atoms_df["type"] = atoms_df["type"].astype(int)

            yield timestep, box_bounds, atoms_df


def calibrate_baseline(dump_path, reference_frame_index=0, threshold_std_multiplier=4.0):
    """
    Calibrates an ABSOLUTE von Mises threshold from a reference frame
    (default: the first frame, i.e. before real loading has produced
    any genuine plastic localization) rather than a per-frame
    percentile.

    THIS REPLACES PERCENTILE-BASED THRESHOLDING, AFTER IT FAILED ON
    REAL DATA. First real-data test (composition 5) reported 15-31
    "bands" at EVERY frame, including timestep=20000 (strain=0.000,
    before any loading at all) -- a percentile threshold always
    returns some "hot" fraction regardless of whether real shear
    localization exists, because it's relative to whatever the
    current frame's distribution looks like, not to an undeformed
    baseline. It was measuring elastic-regime stress fluctuation as
    if it were shear banding. The synthetic self-test never caught
    this because it used near-zero background noise everywhere
    except the engineered band -- unrealistic. Real MD data has
    substantial baseline stress fluctuation at every atom, at every
    frame, which the percentile method happily mistook for bands.

    Returns the absolute von Mises threshold: mean + k*std of the
    reference frame's distribution, so later frames need to show
    stress genuinely exceeding the undeformed baseline, not merely
    rank in the local top percentile of whatever's happening now.
    """
    frames = parse_dump(dump_path)
    for i, (timestep, box_bounds, atoms_df) in enumerate(frames):
        if i == reference_frame_index:
            vm = von_mises_stress(
                atoms_df["sxx"], atoms_df["syy"], atoms_df["szz"],
                atoms_df["sxy"], atoms_df["sxz"], atoms_df["syz"],
            )
            threshold = vm.mean() + threshold_std_multiplier * vm.std()
            print(f"Calibrated from frame {i} (timestep {timestep}): "
                  f"mean={vm.mean():.3e}, std={vm.std():.3e}, threshold={threshold:.3e}")
            return threshold
    raise ValueError(f"reference_frame_index={reference_frame_index} not found in dump")


def detect_bands_in_frame(atoms_df, threshold=None, percentile_threshold=None, eps=6.0, min_samples=15):
    """
    Identify shear bands in a single frame via DBSCAN clustering of
    high-von-Mises atoms.

    Pass an ABSOLUTE threshold (from calibrate_baseline) for real
    data -- this is now the recommended path, not percentile_threshold.
    percentile_threshold is kept only for the synthetic self-test
    below, where the engineered signal-to-noise ratio is controlled
    and known; it should NOT be used on real data given what actually
    happened when it was tried (see calibrate_baseline's docstring).
    """
    vm = von_mises_stress(
        atoms_df["sxx"], atoms_df["syy"], atoms_df["szz"],
        atoms_df["sxy"], atoms_df["sxz"], atoms_df["syz"],
    )
    if threshold is not None:
        hot_mask = vm >= threshold
    elif percentile_threshold is not None:
        hot_mask = vm >= np.percentile(vm, percentile_threshold)
    else:
        raise ValueError("Must supply either threshold (absolute) or percentile_threshold")

    hot_atoms = atoms_df[hot_mask]

    if len(hot_atoms) < min_samples:
        return [], vm

    positions = hot_atoms[["x", "y", "z"]].values
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(positions)

    bands = []
    for label in set(labels):
        if label == -1:  # DBSCAN noise, not a band
            continue
        cluster_mask = labels == label
        cluster_pos = positions[cluster_mask]
        cluster_vm = vm[hot_mask].values[cluster_mask]  # von Mises of this cluster's own atoms
        centroid = cluster_pos.mean(axis=0)
        # PCA-based extent: largest eigenvalue direction approximates
        # band length; spread perpendicular to it approximates width.
        centered = cluster_pos - centroid
        cov = np.cov(centered.T)
        eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
        bands.append({
            "n_atoms": len(cluster_pos),
            "mean_vm": float(cluster_vm.mean()),
            "centroid": centroid,
            "extent_primary": float(np.sqrt(eigvals[0])) if len(eigvals) > 0 else 0.0,
            "extent_secondary": float(np.sqrt(eigvals[1])) if len(eigvals) > 1 else 0.0,
        })

    return bands, vm


def analyze_trajectory(dump_path, eps=6.0, min_samples=15,
                        arrest_growth_tolerance=0.05, threshold_std_multiplier=4.0,
                        erate_per_ps=ERATE_PER_PS):
    """
    Calibrates an absolute von Mises threshold from the first frame
    (the undeformed/elastic baseline), then runs detect_bands_in_frame
    across every frame using that SAME fixed threshold -- not a fresh
    percentile per frame. This is the corrected approach; see
    calibrate_baseline's docstring for why the percentile version
    failed on real data.

    Tracks the largest band's mean von Mises stress vs. strain, with
    the same arrest heuristic as before (plateau detection) -- that
    part was already validated correctly on synthetic data and isn't
    what broke on real data.

    erate_per_ps: MUST be set explicitly for any trajectory that
    doesn't use the baseline 1e9 s^-1 strain rate (erate=0.001 /ps).
    Left as a hardcoded default for a while, which silently produced
    a strain column inflated by exactly the ratio of true rate to
    baseline rate on the first real non-baseline-rate file this was
    run against (rate1e8_1: true final strain 0.15, reported 1.50,
    a 10x error matching 1e9/1e8 exactly) -- band detection itself
    was unaffected (it never used the strain label), but the
    reported strain axis was wrong. Pass the actual rate here for
    any file that isn't a baseline (1e9 /s) run, e.g.
    erate_per_ps=0.0001 for a 1e8 /s trajectory (erate_per_ps =
    strain_rate_per_s / 1e12).
    """
    threshold = calibrate_baseline(dump_path, reference_frame_index=0,
                                    threshold_std_multiplier=threshold_std_multiplier)

    records = []
    for timestep, box_bounds, atoms_df in parse_dump(dump_path):
        bands, vm = detect_bands_in_frame(atoms_df, threshold=threshold, eps=eps, min_samples=min_samples)
        strain = frame_strain(timestep, erate_per_ps=erate_per_ps)
        largest = max(bands, key=lambda b: b["n_atoms"]) if bands else None
        records.append({
            "timestep": timestep,
            "strain": strain,
            "n_bands": len(bands),
            "largest_band_atoms": largest["n_atoms"] if largest else 0,
            "largest_band_mean_vm": largest["mean_vm"] if largest else 0.0,
            "largest_band_extent_primary": largest["extent_primary"] if largest else 0.0,
        })

    df = pd.DataFrame(records).sort_values("strain").reset_index(drop=True)

    arrest_strain = None
    for i in range(1, len(df)):
        prev, curr = df.loc[i - 1, "largest_band_mean_vm"], df.loc[i, "largest_band_mean_vm"]
        if prev > 0 and curr > 0:
            growth = (curr - prev) / prev
            if growth <= arrest_growth_tolerance:
                arrest_strain = df.loc[i, "strain"]
                break

    return df, arrest_strain


# ============================================================
# Self-test against synthetic data built to mimic the real dump
# format exactly, with an engineered high-stress cluster that
# grows then plateaus -- a proxy for shear-band nucleation and
# arrest, so the pipeline's logic can be verified before pointing
# it at real data.
# ============================================================
def _write_synthetic_dump(path, n_atoms=4000, n_frames=10, box_len=None, seed=0):
    """
    box_len defaults to matching the real system's atomic density
    (~0.044 atoms/Ang^3, from 389344 atoms in a ~207 Ang box) rather
    than an arbitrary value.

    IMPORTANT: every atom gets realistic baseline stress noise at
    EVERY frame (not just near-zero background outside the band) --
    the first version of this test didn't, which is exactly why it
    couldn't catch the percentile-threshold bug that real data
    exposed (27 "bands" detected at strain=0.000). A test with
    unrealistically clean background noise can't validate a
    threshold method's ability to reject noise, because there's
    nothing to reject.
    """
    if box_len is None:
        real_density = 389344 / (207.0 ** 3)  # atoms / Ang^3
        box_len = (n_atoms / real_density) ** (1.0 / 3.0)

    rng = np.random.default_rng(seed)
    base_positions = rng.uniform(0, box_len, size=(n_atoms, 3))

    band_slab_frac = 0.15
    band_idx = np.where(base_positions[:, 0] < band_slab_frac * box_len)[0]
    if len(band_idx) < 20:
        band_idx = rng.choice(n_atoms, size=40, replace=False)

    with open(path, "w") as f:
        for frame in range(n_frames):
            timestep = T0_STEPS + frame * 5000
            # Realistic elastic-regime baseline: substantial, present
            # at every frame including frame 0, not near-zero.
            baseline_stress = rng.normal(0, 5.0, size=(n_atoms, 6))

            ramp = min(frame, n_frames // 2) / (n_frames // 2)
            band_stress_magnitude = 15.0 * ramp
            band_stress = rng.normal(band_stress_magnitude, 1.0, size=(len(band_idx), 6))

            stress = baseline_stress.copy()
            stress[band_idx] = band_stress

            f.write("ITEM: TIMESTEP\n")
            f.write(f"{timestep}\n")
            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(f"{n_atoms}\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"0.0 {box_len}\n0.0 {box_len}\n0.0 {box_len}\n")
            f.write("ITEM: ATOMS id type x y z sxx syy szz sxy sxz syz\n")
            for i in range(n_atoms):
                x, y, z = base_positions[i]
                sxx, syy, szz, sxy, sxz, syz = stress[i]
                f.write(f"{i+1} 1 {x:.4f} {y:.4f} {z:.4f} "
                        f"{sxx:.4f} {syy:.4f} {szz:.4f} {sxy:.4f} {sxz:.4f} {syz:.4f}\n")


def _self_test():
    print("=" * 60)
    print("SELF-TEST against synthetic dump -- not real results")
    print("=" * 60)

    synth_path = "synthetic_tension_trajectory.dump"
    _write_synthetic_dump(synth_path, n_atoms=4000, n_frames=10)
    print(f"Wrote synthetic dump: {synth_path}")

    df, arrest_strain = analyze_trajectory(synth_path, eps=8.0, min_samples=10)

    print("\nPer-frame band statistics:")
    print(df.to_string(index=False))

    print(f"\nDetected arrest strain: {arrest_strain}")

    # The critical check this test is actually FOR: zero (or very
    # few) bands at strain=0, since nothing has been engineered
    # there yet -- this is exactly what the old percentile version
    # failed on real data (27 bands reported at strain=0.000).
    n_bands_at_zero_strain = df.loc[df["strain"] == 0.0, "n_bands"].iloc[0]
    print(f"\nBands detected at strain=0.000 (should be 0, or very close): {n_bands_at_zero_strain}")

    if arrest_strain is not None:
        print("PASS: arrest detected, as expected given the engineered plateau.")
    else:
        print("WARNING: no arrest detected -- check DBSCAN parameters "
              "(eps/min_samples/threshold_std_multiplier) against this synthetic case.")

    assert n_bands_at_zero_strain == 0, (
        f"Expected 0 bands at strain=0.000 (nothing engineered there yet), "
        f"got {n_bands_at_zero_strain}. This is exactly the failure mode found "
        f"on real data -- baseline noise is being mistaken for shear bands."
    )
    assert df["n_bands"].max() >= 1, "Expected at least one band once the engineered signal appears."
    detected_phase = df[df["largest_band_mean_vm"] > 0]
    assert len(detected_phase) > 0, "Expected the band to be detected in at least one frame."
    assert detected_phase["largest_band_mean_vm"].std() < detected_phase["largest_band_mean_vm"].mean() * 0.05, (
        "Expected mean von Mises stress to be stable (low relative std) across "
        "frames where the band is actually detected, matching the engineered plateau."
    )
    print("\nAll assertions passed.")


if __name__ == "__main__":
    _self_test()

# ============================================================
# USAGE ON REAL DATA (updated after the percentile-threshold
# approach failed on composition 5's real dump -- 15-31 "bands"
# reported at every frame including strain=0.000, before any
# loading had occurred):
#
#   from detect_shear_bands import analyze_trajectory
#   df, arrest = analyze_trajectory("tension_trajectory.dump")
#
# analyze_trajectory now calibrates an absolute threshold from
# frame 0 (assumed undeformed/elastic) automatically -- no
# percentile_threshold argument anymore. eps/min_samples may still
# need retuning against real cluster sizes (real data's first
# result showed clusters of 20-38 atoms even under the old, too-
# permissive threshold; the new stricter threshold will likely
# produce smaller, cleaner clusters -- check whether min_samples=15
# is still appropriate once real numbers come back, don't assume).
# threshold_std_multiplier (default 4.0) is the other knob to
# revisit: how many standard deviations above the frame-0 baseline
# counts as "genuinely anomalous" is a real modeling choice, not a
# universal constant.
# ============================================================
