"""
interphase_shell_analysis.py

Tests the "weak interphase" hypothesis from external review comment
#7: does yielding/rearrangement nucleate preferentially at the
amorphous/crystalline boundary, rather than being uniformly
distributed through the matrix?

For a given seed center and radius, bins every atom's D2_min
(Falk-Langer non-affine displacement, compute_d2min.py) by its radial
distance from the seed SURFACE (not the center), into:
  - inside the crystal   (r < -shell_width)
  - the interphase shell (-shell_width <= r <= +shell_width), i.e.
    R_seed +/- shell_width, matching the guide's specification
  - the bulk matrix      (r > +shell_width)
and reports the mean D2_min in each region, so an elevated interphase
value (relative to the bulk matrix at the same strain) would support
the weak-interphase hypothesis directly.

THIS SCRIPT DOES NOT FIND SEED CENTERS ITSELF. Locating the crystal
region reliably (so that its true center and radius are known, not
assumed) needs proper structural identification -- Common Neighbor
Analysis or Polyhedral Template Matching, both well-established and
already implemented robustly in OVITO, which is not reproduced here
since a from-scratch reimplementation would be a bigger, more
error-prone undertaking than this binning logic itself, and OVITO is
very likely already part of your workflow given how these dump files
have been handled elsewhere in this project. The recommended
procedure:
  1. Open the composite's FIRST tension-trajectory frame (strain=0,
     before any loading) in OVITO.
  2. Apply Polyhedral Template Matching (or Common Neighbor Analysis)
     to flag HCP-coordinated atoms -- this identifies the inserted
     seed(s) directly, as the only structurally-ordered region in an
     otherwise amorphous matrix.
  3. Select the flagged HCP cluster and read its centroid (Compute
     Property or the cluster analysis modifier's reported center of
     mass) and extent -- this gives seed_center and seed_radius below.
  4. For multi-seed (volume-fraction) composites, repeat per cluster;
     pass each center/radius pair to this script separately, or
     extend the loop below to handle a list of centers at once.

Usage:
    from interphase_shell_analysis import analyze_shell
    result = analyze_shell(dump_path, frame_strain_target=0.071,
                            seed_center=(103.2, 98.7, 101.5), seed_radius=25.0,
                            shell_width=3.0, erate_per_ps=0.001)
"""

import numpy as np
from compute_d2min import parse_dump_frames, compute_d2min_frame, _minimum_image

T0_STEPS = 20000
TIMESTEP_PS = 0.001


def frame_strain(timestep, erate_per_ps):
    return (timestep - T0_STEPS) * erate_per_ps * TIMESTEP_PS


def radial_distance_from_surface(pos, box, seed_center, seed_radius):
    """Signed distance from the seed's SURFACE (not center), using the
    same periodic minimum-image convention as the rest of this
    project's analysis tools -- negative inside the crystal, positive
    in the matrix, zero exactly at the nominal seed surface."""
    L = box[:, 1] - box[:, 0]
    delta = _minimum_image(pos - np.asarray(seed_center), L)
    return np.linalg.norm(delta, axis=1) - seed_radius


def analyze_shell(dump_path, frame_strain_target, seed_center, seed_radius,
                   shell_width=3.0, cutoff=4.5, erate_per_ps=0.001, mode="incremental"):
    """
    mode="cumulative" (the original, and still the default elsewhere in
    this project's docstrings before this fix): compares against frame
    0 always. CONFIRMED UNRELIABLE for this specific shell comparison at
    real strains: on vf_1 at strain ~0.07, the bulk_matrix reference
    came back with mean D2_min=0.678 -- not a clean near-zero baseline,
    but already substantially elevated by the same generic drift effect
    documented in analyze_d2min_trajectory.py (comparing against an
    increasingly distant fixed reference accumulates displacement from
    ordinary thermal motion alone, not just real rearrangement). This
    contaminates the bulk reference the interphase region is compared
    against, undermining the ratio.

    mode="incremental" (now the default): compares each frame only
    against the frame immediately before it, avoiding this drift
    contamination -- the same fix already applied to
    analyze_d2min_trajectory.py, carried over here since this tool was
    originally built without it.
    """
    frames = parse_dump_frames(dump_path)
    strains = np.array([frame_strain(fr[0], erate_per_ps) for fr in frames])
    i = int(np.argmin(np.abs(strains - frame_strain_target)))
    actual_strain = strains[i]

    if mode == "cumulative":
        ref_ids, ref_pos, ref_box = frames[0][1], frames[0][2], frames[0][3]
    else:
        assert i > 0, "incremental mode needs at least one frame before the target"
        ref_ids, ref_pos, ref_box = frames[i-1][1], frames[i-1][2], frames[i-1][3]
    timestep, ids, pos, box = frames[i]
    common, d2min = compute_d2min_frame(ref_ids, ref_pos, ref_box, ids, pos, box, cutoff)
    valid = ~np.isnan(d2min)

    common_pos = pos[np.searchsorted(ids, common)]
    r = radial_distance_from_surface(common_pos, box, seed_center, seed_radius)

    # IMPORTANT: D2_min is a neighborhood-based quantity -- it depends
    # on how an atom's NEIGHBORS moved, not just the atom itself. An
    # atom binned as "just outside" the shell can still have most of
    # its neighbor sphere actually falling inside the shell if binned
    # right at the boundary, contaminating its score with the shell's
    # behavior even though the atom itself is categorized as
    # crystal/matrix. Confirmed directly on synthetic data: without
    # this gap, "crystal interior" came out ~25x inflated purely from
    # boundary atoms whose neighbor shells reached into an engineered
    # shell-only rearrangement. A gap of one cutoff radius on each
    # side of the shell guarantees no neighbor-sphere overlap; atoms
    # in the gap itself are ambiguous and excluded from all three
    # categories rather than mislabeled.
    regions = {
        "inside_crystal": (r < -(shell_width + cutoff)) & valid,
        "interphase_shell": (r >= -shell_width) & (r <= shell_width) & valid,
        "bulk_matrix": (r > (shell_width + cutoff)) & valid,
    }
    result = dict(target_strain=frame_strain_target, actual_strain=actual_strain, regions={})
    for name, mask in regions.items():
        n = int(mask.sum())
        result["regions"][name] = dict(
            n_atoms=n,
            mean_d2min=float(np.nanmean(d2min[mask])) if n else float("nan"),
            median_d2min=float(np.nanmedian(d2min[mask])) if n else float("nan"),
        )
    return result


def print_result(result):
    print(f"Strain requested {result['target_strain']:.3f}, actual frame strain {result['actual_strain']:.3f}")
    for name, r in result["regions"].items():
        print(f"  {name:<18} n={r['n_atoms']:>7}  mean D2_min={r['mean_d2min']:.4f}  "
              f"median D2_min={r['median_d2min']:.4f}")
    shell = result["regions"]["interphase_shell"]["mean_d2min"]
    bulk = result["regions"]["bulk_matrix"]["mean_d2min"]
    if bulk > 0:
        print(f"  Interphase/bulk ratio: {shell / bulk:.2f}x "
              f"({'ELEVATED -- consistent with weak-interphase hypothesis' if shell > bulk * 1.2 else 'not clearly elevated'})")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("dump_path")
        parser.add_argument("target_strain", type=float)
        parser.add_argument("center_x", type=float)
        parser.add_argument("center_y", type=float)
        parser.add_argument("center_z", type=float)
        parser.add_argument("--seed_radius", type=float, required=True,
                             help="Found via find_seed_centers.py. NOTE: that tool's radius "
                                  "estimate is a systematic UNDERESTIMATE of the true seed size "
                                  "(confirmed on vf_1: recovered ~9.7 Ang against a known 12 Ang "
                                  "build target, ~20% low) -- surface atoms of the seed partially "
                                  "blend with the matrix and fail the strict candidate threshold, "
                                  "shrinking the apparent cluster. Consider running this with BOTH "
                                  "the raw find_seed_centers.py value and the known nominal build "
                                  "radius (e.g. 12.0 for the baseline VF sweep) to bracket the true "
                                  "interphase location, rather than trusting either alone.")
        parser.add_argument("--shell_width", type=float, default=3.0)
        parser.add_argument("--cutoff", type=float, default=4.5)
        parser.add_argument("--erate_per_ps", type=float, default=0.001)
        parser.add_argument("--mode", choices=["cumulative", "incremental"], default="incremental",
                             help="incremental (default) avoids the drift contamination found in "
                                  "the cumulative mode at real strains; pass --mode cumulative only "
                                  "to reproduce the earlier, now-understood-to-be-unreliable result.")
        args = parser.parse_args()
        result = analyze_shell(args.dump_path, args.target_strain,
                                (args.center_x, args.center_y, args.center_z), args.seed_radius,
                                shell_width=args.shell_width, cutoff=args.cutoff,
                                erate_per_ps=args.erate_per_ps, mode=args.mode)
        print_result(result)
        sys.exit(0)

    print("Self-test: synthetic composite with a known crystal region and an "
          "engineered interphase-only rearrangement.\n")

    rng = np.random.default_rng(11)
    n, L = 6000, 60.0
    origin = np.array([5.0, -10.0, 20.0])
    pos = rng.uniform(0, L, size=(n, 3)) + origin
    seed_center = origin + np.array([30.0, 30.0, 30.0])
    seed_radius = 15.0

    r0 = radial_distance_from_surface(pos, np.array([origin, origin + L, origin + L]).T, seed_center, seed_radius) \
        if False else None  # (box arg shape handled directly below instead)
    box = np.column_stack([origin, origin + L])
    r0 = radial_distance_from_surface(pos, box, seed_center, seed_radius)

    # Engineer a clean, unambiguous test: ONLY atoms in the interphase
    # shell (-3 to +3 Ang from the nominal surface) get a real
    # non-affine kick; atoms deep in the "crystal" and deep in the
    # "matrix" get only thermal noise. Built as a PROPERLY ACCUMULATING
    # trajectory (each frame = previous frame + one small independent
    # step), not independent redraws from the pristine base -- the
    # latter was tried first and confirmed to inflate the incremental
    # comparison's baseline variance by exactly 2x (two independent
    # noise draws compared against each other, rather than one new
    # step building on the last), the same issue already found and
    # fixed once before in this project's sticky-defect D2min test.
    shell_mask = (r0 >= -3.0) & (r0 <= 3.0)
    print(f"Synthetic interphase-shell atom count: {shell_mask.sum()}")
    kick_dir = rng.normal(size=3); kick_dir /= np.linalg.norm(kick_dir)

    frame0 = pos.copy()
    frame1 = frame0 + rng.normal(scale=0.03, size=(n, 3))
    frame2 = frame1 + rng.normal(scale=0.03, size=(n, 3))
    frame2[shell_mask] += kick_dir * 1.0

    ids = np.arange(1, n + 1)
    tmp_dump = "synthetic_interphase.dump"
    with open(tmp_dump, "w") as f:
        for ts, p in [(20000, frame0), (25000, frame1), (30000, frame2)]:
            f.write("ITEM: TIMESTEP\n"); f.write(f"{ts}\n")
            f.write("ITEM: NUMBER OF ATOMS\n"); f.write(f"{n}\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"{origin[0]} {origin[0]+L}\n{origin[1]} {origin[1]+L}\n{origin[2]} {origin[2]+L}\n")
            f.write("ITEM: ATOMS id type x y z\n")
            for i in range(n):
                f.write(f"{ids[i]} 1 {p[i,0]:.5f} {p[i,1]:.5f} {p[i,2]:.5f}\n")

    result = analyze_shell(tmp_dump, frame_strain_target=0.01, seed_center=seed_center,
                            seed_radius=seed_radius, shell_width=3.0, cutoff=4.0, mode="incremental")
    print_result(result)

    shell_d2 = result["regions"]["interphase_shell"]["mean_d2min"]
    bulk_d2 = result["regions"]["bulk_matrix"]["mean_d2min"]
    crystal_d2 = result["regions"]["inside_crystal"]["mean_d2min"]
    print(f"\n  (gap-buffered) bulk mean = {bulk_d2:.5f}, crystal mean = {crystal_d2:.5f} "
          f"-- both should now be near the pure-thermal-noise floor (~0.003), not elevated.")
    assert shell_d2 > 15 * bulk_d2, "FAIL: engineered interphase rearrangement should be sharply elevated vs. bulk"
    assert shell_d2 > 15 * crystal_d2, "FAIL: engineered interphase rearrangement should be sharply elevated vs. crystal interior"
    assert crystal_d2 < 0.01, "FAIL: crystal interior should be near the thermal-noise floor once properly buffered from the shell"
    print("PASS: interphase-only rearrangement correctly isolated from both the crystal interior and the bulk matrix, "
          "with both baselines now near the thermal-noise floor once properly buffered.")

    import os
    os.remove(tmp_dump)
