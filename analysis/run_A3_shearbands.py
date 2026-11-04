#!/usr/bin/env python3
"""run_A3_shearbands.py -- DBSCAN von Mises shear-band check (Section 2.6 settings, unchanged) on the new tension dumps.

  python3 run_A3_shearbands.py glass0_r5 glass0_r1 glass0_r6 vf_r5_1 vf_r5_2 vf_r5_3 temp_350 temp_400 rate2e9_run4
  python3 run_A3_shearbands.py rep_r5_a rep_r5_b ... temp_250 temp_300 rate2e9_run3

Each argument is a case folder containing tension_trajectory.dump.  Settings are the ones used for the original 25 runs
(eps 6.0 A, min_samples 15, threshold = frame-0 mean + 4 sd).  The strain axis uses erate 1e-3/ps, except cases whose
name contains 'rate2e9' (2e-3/ps).  Writes <case>/A3_shear_frames.csv (per frame) and A3_shear_summary.csv (one row per case).
hcp_x / hcp_z are excluded (no amorphous matrix; stress field is that of a crystal)."""
import sys, os, time
import pandas as pd
from detect_shear_bands import analyze_trajectory
rows = []
for case in sys.argv[1:]:
    dump = os.path.join(case, "tension_trajectory.dump")
    if not os.path.exists(dump):
        print("MISSING", dump); continue
    t0 = time.time()
    erate = 0.002 if "rate2e9" in case else 0.001
    df, arrest = analyze_trajectory(dump, erate_per_ps=erate)
    df.to_csv(os.path.join(case, "A3_shear_frames.csv"), index=False)
    imax = df.n_bands.idxmax()
    late = df[df.strain >= 0.10]
    rows.append(dict(case=case, erate_per_ps=erate, frames=len(df), max_bands=int(df.n_bands.max()), strain_at_max_bands=float(df.strain[imax]),
                     bands_at_frame0=int(df.n_bands.iloc[0]), largest_band_atoms=int(df.largest_band_atoms.max()),
                     mean_bands_strain_ge_10pct=float(late.n_bands.mean()) if len(late) else float("nan"), arrest_strain=arrest,
                     seconds=round(time.time() - t0)))
    print(rows[-1], flush=True)
out = pd.DataFrame(rows); out.to_csv("A3_shear_summary.csv", index=False)
print(out.to_string())
