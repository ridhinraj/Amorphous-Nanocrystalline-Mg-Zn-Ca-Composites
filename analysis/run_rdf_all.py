"""
run_rdf_all.py

Runs compute_rdf.py's logic across all 10 composition directories
(not vf_1-3, which share run_1's matrix structure and would just
reproduce run_1's RDF) and writes one consolidated summary table
alongside the full per-composition RDF CSVs, so Figure 2 and the
Section 3.2 comparison table can be built from a single pass instead
of ten separate manual runs.

Run from the project root:
    source ~/shap_env/bin/activate
    python3 run_rdf_all.py
"""

import os
import numpy as np
import pandas as pd

from compute_rdf import read_lammps_data, compute_gr

COMPOSITIONS = {
    "run_1": (0.35, 0.05), "run_2": (0.30, 0.04), "run_3": (0.40, 0.06),
    "run_4": (0.25, 0.05), "run_5": (0.35, 0.03), "run_6": (0.35, 0.07),
    "run_7": (0.30, 0.06), "run_8": (0.40, 0.04), "run_9": (0.32, 0.05),
    "run_10": (0.38, 0.05),
}
ATOMIC_MASSES = {1: 24.305, 2: 65.38, 3: 40.078}  # Mg, Zn, Ca


def main():
    summary_rows = []

    for run_id, (xZn, xCa) in COMPOSITIONS.items():
        data_path = os.path.join(run_id, "quenched_matrix_300K.data")
        if not os.path.exists(data_path):
            print(f"{run_id}: SKIPPED (no quenched_matrix_300K.data found)")
            summary_rows.append(dict(run_id=run_id, xZn=xZn, xCa=xCa, status="missing"))
            continue

        try:
            positions, types, box_lengths = read_lammps_data(data_path)
            box_length = box_lengths[0]
            N = len(positions)

            r, g_total = compute_gr(positions, box_length, r_max=10.0, n_bins=200)

            mask = r > 1.5
            peak_idx = np.argmax(g_total[mask])
            peak_r = r[mask][peak_idx]
            peak_g = g_total[mask][peak_idx]
            g_at_rmax = g_total[-1]

            xMg = 1.0 - xZn - xCa
            m_avg = xMg * ATOMIC_MASSES[1] + xZn * ATOMIC_MASSES[2] + xCa * ATOMIC_MASSES[3]
            density_g_cm3 = (N * m_avg / 6.02214076e23) / (box_length ** 3 * 1e-24)

            out_csv = os.path.join(run_id, "quenched_matrix_300K_rdf.csv")
            pd.DataFrame({"r": r, "g_total": g_total}).to_csv(out_csv, index=False)

            print(f"{run_id}: N={N}, box={box_length:.3f} Ang, "
                  f"peak r={peak_r:.3f} Ang, peak g={peak_g:.3f}, "
                  f"g(r_max)={g_at_rmax:.3f}, density={density_g_cm3:.3f} g/cm^3")

            summary_rows.append(dict(
                run_id=run_id, xZn=xZn, xCa=xCa, N=N, box_length=box_length,
                first_peak_r=peak_r, first_peak_g=peak_g, g_at_rmax=g_at_rmax,
                density_g_cm3=density_g_cm3, status="OK",
            ))
        except Exception as e:
            print(f"{run_id}: ERROR: {e}")
            summary_rows.append(dict(run_id=run_id, xZn=xZn, xCa=xCa, status=f"ERROR: {e}"))

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv("rdf_summary_all_compositions.csv", index=False)
    print(f"\nWritten: rdf_summary_all_compositions.csv ({len(summary_df)} rows)")

    n_ok = (summary_df["status"] == "OK").sum()
    print(f"{n_ok} of {len(summary_df)} compositions processed successfully.")

    if n_ok > 0:
        ok = summary_df[summary_df["status"] == "OK"]
        print(f"\ng(r_max) range: {ok['g_at_rmax'].min():.3f} to {ok['g_at_rmax'].max():.3f} "
              f"(all should be close to 1.0 for confirmed amorphous character)")
        print(f"First-peak position range: {ok['first_peak_r'].min():.3f} to {ok['first_peak_r'].max():.3f} Ang")
        print(f"Density range: {ok['density_g_cm3'].min():.3f} to {ok['density_g_cm3'].max():.3f} g/cm^3")


if __name__ == "__main__":
    main()
