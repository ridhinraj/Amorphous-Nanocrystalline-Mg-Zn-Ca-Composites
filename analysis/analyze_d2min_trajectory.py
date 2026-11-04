"""
analyze_d2min_trajectory.py

Runs D2_min (compute_d2min.py) across a full tension trajectory,
comparing every recorded frame against the UNDEFORMED FIRST FRAME --
the same reference-frame convention already used by the von Mises
threshold calibration in detect_shear_bands.py, so the two analyses
are directly comparable frame-for-frame.

Reports, per frame: the strain, the number of atoms with D2_min above
an absolute threshold (calibrated the same way as the von Mises
threshold: from the natural spread of D2_min values in the reference
frame's own thermal noise, via a short probe comparison of frame 0
against frame 1), and whether those high-D2_min atoms cluster
spatially (DBSCAN, matching detect_shear_bands.py's clustering step).

This is a genuine cross-check, not a re-labelling of the same result:
von Mises stress and D2_min respond to different things (elastic
stress vs. irreversible rearrangement), so agreement between the two
is real evidence, not a foregone conclusion.

Usage:
    python3 analyze_d2min_trajectory.py <dump_path> [--erate_per_ps 0.001]
"""

import sys
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from compute_d2min import parse_dump_frames, compute_d2min_frame

T0_STEPS = 20000
TIMESTEP_PS = 0.001


def frame_strain(timestep, erate_per_ps):
    return (timestep - T0_STEPS) * erate_per_ps * TIMESTEP_PS


def calibrate_d2min_threshold(frames, cutoff, std_multiplier=4.0):
    """Same philosophy as detect_shear_bands.py's calibrate_baseline:
    probe the reference frame's own short-timescale thermal noise
    (frame 0 vs frame 1, i.e. 5000 steps of pure elastic/thermal
    motion before real plastic strain has accumulated) rather than
    assuming a fixed threshold, since the right scale depends on the
    system."""
    ref_ids, ref_pos, ref_box = frames[0][1], frames[0][2], frames[0][3]
    probe_ids, probe_pos, probe_box = frames[1][1], frames[1][2], frames[1][3]
    _, d2min_probe = compute_d2min_frame(ref_ids, ref_pos, ref_box, probe_ids, probe_pos, probe_box, cutoff)
    valid = d2min_probe[~np.isnan(d2min_probe)]
    mean, std = valid.mean(), valid.std()
    threshold = mean + std_multiplier * std
    print(f"D2_min threshold calibrated from frame 0->1 (thermal noise only): "
          f"mean={mean:.4f}, std={std:.4f}, threshold={threshold:.4f}")
    return threshold


def analyze(dump_path, erate_per_ps=0.001, cutoff=4.5, eps=6.0, min_samples=15, mode="cumulative"):
    """
    mode="cumulative": every frame compared against the UNDEFORMED
    FIRST FRAME (frame 0), matching the convention already used for
    the von Mises threshold in detect_shear_bands.py. An atom flagged
    at some strain typically STAYS flagged at every later strain too,
    since displacement accumulated from a fixed reference can only
    grow -- this is a real property of the measurement, not a bug,
    but it means late-strain "coalescence into one big cluster" could
    partly reflect this bookkeeping rather than deformation actually
    becoming less localized in real time.

    mode="incremental": each frame compared only against the frame
    immediately BEFORE it (a single 5000-step interval throughout),
    which measures what is rearranging RIGHT NOW at each strain step
    and does not carry forward atoms that rearranged earlier and then
    stopped. This is the version that can distinguish "these atoms
    have EVER rearranged" from "these atoms are CURRENTLY rearranging".
    Both modes use the same threshold, calibrated once from the
    frame-0-to-1 thermal-noise-only interval, since every incremental
    comparison and the calibration interval itself span the same
    5000 steps.
    """
    frames = parse_dump_frames(dump_path)
    threshold = calibrate_d2min_threshold(frames, cutoff)

    ref0_ids, ref0_pos, ref0_box = frames[0][1], frames[0][2], frames[0][3]

    records = []
    prev = frames[0]
    for k, (timestep, ids, pos, box) in enumerate(frames):
        if mode == "cumulative":
            r_ids, r_pos, r_box = ref0_ids, ref0_pos, ref0_box
        else:
            r_ids, r_pos, r_box = prev[1], prev[2], prev[3]

        common, d2min = compute_d2min_frame(r_ids, r_pos, r_box, ids, pos, box, cutoff)
        valid = ~np.isnan(d2min)
        high = valid & (d2min > threshold)
        n_bands, largest = 0, 0
        if high.sum() >= min_samples:
            coords = pos[np.searchsorted(ids, common[high])]
            labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords)
            real = labels[labels >= 0]
            if len(real):
                n_bands = len(set(real))
                largest = int(np.bincount(real).max())
        records.append(dict(
            timestep=timestep, strain=frame_strain(timestep, erate_per_ps), mode=mode,
            n_high_d2min_atoms=int(high.sum()), n_d2min_clusters=n_bands,
            largest_d2min_cluster=largest,
            mean_d2min=float(np.nanmean(d2min[valid])),
        ))
        print(f"  frame {timestep}: strain={records[-1]['strain']:.3f}, "
              f"n_high={records[-1]['n_high_d2min_atoms']}, clusters={n_bands}, "
              f"largest={largest}")
        prev = (timestep, ids, pos, box)

    return pd.DataFrame(records)


if __name__ == "__main__":
    dump_path = sys.argv[1] if len(sys.argv) > 1 else "tension_trajectory.dump"
    erate = 0.001
    if "--erate_per_ps" in sys.argv:
        erate = float(sys.argv[sys.argv.index("--erate_per_ps") + 1])
    mode = "cumulative"
    if "--mode" in sys.argv:
        mode = sys.argv[sys.argv.index("--mode") + 1]
        assert mode in ("cumulative", "incremental"), "mode must be 'cumulative' or 'incremental'"
    df = analyze(dump_path, erate_per_ps=erate, mode=mode)
    out = dump_path.replace(".dump", f"_d2min_{mode}_summary.csv").replace("/", "__")
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")
    print(f"Max clusters at any frame: {df['n_d2min_clusters'].max()}")
