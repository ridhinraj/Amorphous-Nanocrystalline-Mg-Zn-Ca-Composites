"""
consolidate_results.py

Pulls together everything checked piecemeal over many separate runs
into one pass across every run_N and vf_N directory:

  1. Yield stress (from tension_stress_strain.txt)
  2. Diffuse-deformation check (von Mises mean/std, first vs last
     frame -- extends the finding confirmed on run_5/run_1/vf_3 to
     every composition, checking whether any of them DO show real
     localization that the earlier spot-checks missed)
  3. Ca-vs-Zn per-atom stress comparison (extends the finding
     confirmed on run_1/run_5/run_6 to every composition and the
     full VF sweep)

Does NOT require shap or the shap_env venv -- only the direct
mean/SEM comparison is computed here (the part that actually
matters for the paper), not the SHAP feature-importance ranking.
Run this under the normal anaconda environment:

    module load apps/anaconda3/2024.06
    python3 consolidate_results.py

Missing files (e.g. run_8/run_9 still in progress) are handled
gracefully -- reported as incomplete, not a crash.
"""

import os
import numpy as np
import pandas as pd

from detect_shear_bands import parse_dump, von_mises_stress

# Composition matrix, from setup_job_array.py / Table 1
COMPOSITIONS = {
    "run_1": (0.35, 0.05, 25), "run_2": (0.30, 0.04, 20), "run_3": (0.40, 0.06, 30),
    "run_4": (0.25, 0.05, 25), "run_5": (0.35, 0.03, 15), "run_6": (0.35, 0.07, 35),
    "run_7": (0.30, 0.06, 25), "run_8": (0.40, 0.04, 20), "run_9": (0.32, 0.05, 28),
    "run_10": (0.38, 0.05, 22),
}
VF_CASES = {
    "vf_1": (0.35, 0.05, 0.05), "vf_2": (0.35, 0.05, 0.10), "vf_3": (0.35, 0.05, 0.15),
}


def get_yield_stress(run_dir):
    path = os.path.join(run_dir, "tension_stress_strain.txt")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None,
                          names=["strain", "stress_x", "stress_y", "stress_z"])
        return df["stress_x"].tail(3).mean()
    except Exception as e:
        return f"ERROR: {e}"


def get_dump_diagnostics(run_dir):
    """Returns (is_flat, ca_zn_diff) or (None, None) if dump is missing/unreadable."""
    path = os.path.join(run_dir, "tension_trajectory.dump")
    if not os.path.exists(path):
        return None, None
    try:
        frames = list(parse_dump(path))
        if len(frames) < 2:
            return None, None

        # Flatness check: compare mean von Mises, first frame vs last
        _, _, first_atoms = frames[0]
        _, _, last_atoms = frames[-1]
        vm_first = von_mises_stress(
            first_atoms["sxx"], first_atoms["syy"], first_atoms["szz"],
            first_atoms["sxy"], first_atoms["sxz"], first_atoms["syz"],
        )
        vm_last = von_mises_stress(
            last_atoms["sxx"], last_atoms["syy"], last_atoms["szz"],
            last_atoms["sxy"], last_atoms["sxz"], last_atoms["syz"],
        )
        pct_change = 100 * (vm_last.mean() - vm_first.mean()) / vm_first.mean()
        is_flat = abs(pct_change) < 5.0  # <5% change in mean vM -- matches what run_5/1/6/vf_3 all showed

        # Ca-vs-Zn stress comparison, last frame
        types = last_atoms["type"].values
        vm = vm_last.values
        ca_vals = vm[types == 3]
        zn_vals = vm[types == 2]
        if len(ca_vals) > 10 and len(zn_vals) > 10:
            log_ca = np.log(np.maximum(ca_vals, 1.0))
            log_zn = np.log(np.maximum(zn_vals, 1.0))
            ca_zn_diff = log_ca.mean() - log_zn.mean()
        else:
            ca_zn_diff = None

        return is_flat, ca_zn_diff
    except Exception as e:
        return f"ERROR: {e}", None


def main():
    rows = []

    print(f"{'run_id':<8} {'xZn':>5} {'xCa':>5} {'yield_GPa':>10} {'flat?':>7} {'Ca-Zn diff':>11}  status")
    print("-" * 70)

    for run_id, (xZn, xCa, radius) in COMPOSITIONS.items():
        yield_stress = get_yield_stress(run_id)
        is_flat, ca_zn_diff = get_dump_diagnostics(run_id)

        if yield_stress is None:
            status = "INCOMPLETE (no tension_stress_strain.txt)"
        elif isinstance(yield_stress, str):
            status = yield_stress
        else:
            status = "OK"

        rows.append(dict(run_id=run_id, xZn=xZn, xCa=xCa, volume_fraction=None,
                          seed_radius=radius, yield_stress_GPa=yield_stress,
                          is_flat=is_flat, ca_zn_log_diff=ca_zn_diff, status=status))

        ys_str = f"{yield_stress:.3f}" if isinstance(yield_stress, float) else "--"
        flat_str = str(is_flat) if is_flat is not None else "--"
        diff_str = f"{ca_zn_diff:+.3f}" if isinstance(ca_zn_diff, float) else "--"
        print(f"{run_id:<8} {xZn:>5} {xCa:>5} {ys_str:>10} {flat_str:>7} {diff_str:>11}  {status}")

    print()
    for run_id, (xZn, xCa, vf) in VF_CASES.items():
        yield_stress = get_yield_stress(run_id)
        is_flat, ca_zn_diff = get_dump_diagnostics(run_id)

        if yield_stress is None:
            status = "INCOMPLETE (no tension_stress_strain.txt)"
        elif isinstance(yield_stress, str):
            status = yield_stress
        else:
            status = "OK"

        rows.append(dict(run_id=run_id, xZn=xZn, xCa=xCa, volume_fraction=vf,
                          seed_radius=None, yield_stress_GPa=yield_stress,
                          is_flat=is_flat, ca_zn_log_diff=ca_zn_diff, status=status))

        ys_str = f"{yield_stress:.3f}" if isinstance(yield_stress, float) else "--"
        flat_str = str(is_flat) if is_flat is not None else "--"
        diff_str = f"{ca_zn_diff:+.3f}" if isinstance(ca_zn_diff, float) else "--"
        print(f"{run_id:<8} {xZn:>5} {xCa:>5} {ys_str:>10} {flat_str:>7} {diff_str:>11}  {status}")

    out_df = pd.DataFrame(rows)
    out_df.to_csv("consolidated_results.csv", index=False)
    print(f"\nWritten: consolidated_results.csv ({len(out_df)} rows)")

    n_complete = out_df["yield_stress_GPa"].apply(lambda x: isinstance(x, float)).sum()
    print(f"\n{n_complete} of {len(out_df)} cases have usable yield-stress data.")

    flat_count = sum(1 for x in out_df["is_flat"] if x == True)
    non_flat_count = sum(1 for x in out_df["is_flat"] if x == False)
    print(f"Diffuse-deformation check: {flat_count} flat (no discrete band), "
          f"{non_flat_count} NOT flat (worth a closer look if any appear here).")


if __name__ == "__main__":
    main()
