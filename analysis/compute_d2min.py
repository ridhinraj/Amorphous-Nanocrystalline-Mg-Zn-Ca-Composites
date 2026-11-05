"""
compute_d2min.py

Falk & Langer (1998) non-affine squared displacement D2_min, computed
between a REFERENCE frame (typically the undeformed, zero-strain
frame) and a CURRENT frame, for every atom in a LAMMPS dump
trajectory.

WHY THIS IS NEEDED (per external review): per-atom von Mises stress,
used throughout this project's deformation-character analysis
(detect_shear_bands.py) and solute-clustering analysis
(solute_clustering_shap.py), cannot distinguish reversible ELASTIC
stress concentration from irreversible PLASTIC rearrangement -- an
atom can carry high stress without ever having actually rearranged,
and vice versa. D2_min is the standard, purpose-built measure for
irreversible local rearrangement: for each atom, it fits the best
local affine deformation (a 3x3 matrix F) that maps its neighbors'
reference-frame relative positions onto their current-frame relative
positions, and reports the mean-squared RESIDUAL after removing that
best-fit affine part. An atom whose whole neighborhood moved by a
single affine transformation (i.e. just elastic strain) has D2_min
near zero; an atom that rearranged non-affinely (a plastic event) has
a large D2_min, regardless of how much stress it happens to carry at
that instant.

Reference: M.L. Falk & J.S. Langer, "Dynamics of viscoplastic
deformation in amorphous solids," Phys. Rev. E 57, 7192 (1998).

Usage:
    from compute_d2min import parse_dump_frames, compute_d2min_frame

    frames = parse_dump_frames("tension_trajectory.dump")
    ref_ids, ref_pos, ref_box = frames[0]
    for timestep, ids, pos, box in frames[1:]:
        d2min = compute_d2min_frame(ref_ids, ref_pos, ref_box, ids, pos, box, cutoff=4.5)
        # d2min is an array aligned to `ids`, one value per atom
"""

import numpy as np


def parse_dump_frames(dump_path):
    """Parses a LAMMPS custom dump file into a list of
    (timestep, atom_ids, positions[N,3], box_bounds[3,2]) tuples,
    reading id, x, y, z columns (any additional columns are ignored)."""
    frames = []
    with open(dump_path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(f.readline().strip())
            f.readline()  # ITEM: NUMBER OF ATOMS
            n_atoms = int(f.readline().strip())
            f.readline()  # ITEM: BOX BOUNDS ...
            box = np.array([list(map(float, f.readline().split())) for _ in range(3)])
            header = f.readline().strip()  # ITEM: ATOMS id type x y z ...
            cols = header.replace("ITEM: ATOMS", "").split()
            id_idx = cols.index("id")
            xi, yi, zi = cols.index("x"), cols.index("y"), cols.index("z")
            ids = np.empty(n_atoms, dtype=np.int64)
            pos = np.empty((n_atoms, 3), dtype=np.float64)
            for i in range(n_atoms):
                parts = f.readline().split()
                ids[i] = int(parts[id_idx])
                pos[i] = [float(parts[xi]), float(parts[yi]), float(parts[zi])]
            order = np.argsort(ids)
            frames.append((timestep, ids[order], pos[order], box))
    return frames


def _minimum_image(delta, box_lengths):
    """Wraps displacement vectors into [-L/2, L/2) per periodic dimension."""
    return delta - box_lengths * np.round(delta / box_lengths)


def compute_d2min_frame(ref_ids, ref_pos, ref_box, cur_ids, cur_pos, cur_box, cutoff):
    """D2_min for every atom present in both frames, using the
    REFERENCE frame's neighbor list (the standard Falk-Langer choice:
    neighbors are defined once, at the reference configuration, not
    re-found at every frame). Returns an array aligned to a sorted,
    common atom-id ordering.

    ref_box/cur_box: the 3x2 box-bounds array as read by
    parse_dump_frames (used only for minimum-image wrapping; assumes
    an orthogonal box, matching every dump this project has produced).
    """
    common, ref_idx, cur_idx = np.intersect1d(ref_ids, cur_ids, return_indices=True)
    n = len(common)
    rp = ref_pos[ref_idx]
    cp = cur_pos[cur_idx]
    ref_L = ref_box[:, 1] - ref_box[:, 0]
    cur_L = cur_box[:, 1] - cur_box[:, 0]

    # cKDTree's periodic mode requires coordinates in [0, boxsize) --
    # real LAMMPS dumps almost never start their box at the origin
    # (confirmed directly: run_1 failed with "input data greater than
    # the size of the periodic box" on first real use, the same
    # nonzero-origin issue found and fixed in compute_rdf.py earlier
    # in this project). Shift into [0, L) by the box's own lower
    # bound and wrap, rather than assuming an origin-anchored box.
    rp_wrapped = (rp - ref_box[:, 0]) % ref_L

    from scipy.spatial import cKDTree
    tree = cKDTree(rp_wrapped, boxsize=ref_L)

    d2min = np.full(n, np.nan)
    for i in range(n):
        neigh = tree.query_ball_point(rp_wrapped[i], cutoff)
        neigh = [j for j in neigh if j != i]
        if len(neigh) < 4:
            continue  # too few neighbors for a stable 3x3 fit
        dref = _minimum_image(rp[neigh] - rp[i], ref_L)
        dcur = _minimum_image(cp[neigh] - cp[i], cur_L)

        X = dcur.T @ dref
        Y = dref.T @ dref
        try:
            F = X @ np.linalg.inv(Y)
        except np.linalg.LinAlgError:
            continue

        predicted = dref @ F.T
        residual = dcur - predicted
        d2min[i] = np.mean(np.sum(residual ** 2, axis=1))

    return common, d2min


if __name__ == "__main__":
    print("Self-test: synthetic affine shear vs. a single non-affine rearrangement.\n")
    rng = np.random.default_rng(0)

    n_atoms = 2000
    L = 40.0
    ref_pos = rng.uniform(0, L, size=(n_atoms, 3))
    ref_box = np.array([[0, L], [0, L], [0, L]])
    ids = np.arange(1, n_atoms + 1)

    # Test 1: pure affine UNIAXIAL stretch applied to every atom --
    # matching the actual use case (tension trajectories only stretch/
    # contract box lengths, never shear/tilt the box). D2_min should
    # be ~0 everywhere, since a perfectly affine map has zero residual
    # by construction, AND periodic wrapping stays valid (an
    # orthogonal box stretched into another orthogonal box), unlike a
    # sheared box, which would need triclinic (tilted) periodic
    # images that this minimum-image implementation does not provide
    # -- deliberately not tested here, since it is not needed for any
    # trajectory this project has produced.
    strain, nu = 0.05, 0.3
    F_affine = np.diag([1 + strain, 1 - nu * strain, 1 - nu * strain])
    cur_pos_affine = ref_pos @ F_affine.T
    cur_box_affine = ref_box * np.array([[1, 1 + strain], [1, 1 - nu * strain], [1, 1 - nu * strain]])
    # (the above scales only the upper bound since ref_box's lower bound is 0; for a
    # box not starting at the origin this would need both bounds scaled from the center)

    _, d2_affine = compute_d2min_frame(ids, ref_pos, ref_box, ids, cur_pos_affine, cur_box_affine, cutoff=4.0)
    print(f"Pure affine uniaxial stretch: mean D2_min = {np.nanmean(d2_affine):.2e}, "
          f"max D2_min = {np.nanmax(d2_affine):.2e}  (both should be ~0)")
    assert np.nanmax(d2_affine) < 1e-20, "FAIL: affine deformation should give exactly zero D2_min"
    print("PASS: affine-only deformation gives D2_min consistent with zero.\n")

    # Test 2: same affine shear, but ONE atom additionally gets a large
    # random (non-affine) kick. D2_min should now be clearly elevated
    # at and near that one atom, and still ~0 far away.
    cur_pos_defect = cur_pos_affine.copy()
    defect_atom = 500
    kick = rng.normal(scale=2.0, size=3)
    # Keep the kicked atom safely away from any periodic boundary in
    # this synthetic test, so the check below isolates the algorithm's
    # response to a genuine local rearrangement rather than mixing in
    # any boundary artifact.
    ref_pos[defect_atom] = np.clip(ref_pos[defect_atom], 8.0, L - 8.0)
    cur_pos_defect[defect_atom] = ref_pos[defect_atom] @ F_affine.T + kick

    _, d2_defect = compute_d2min_frame(ids, ref_pos, ref_box, ids, cur_pos_defect, cur_box_affine, cutoff=4.0)

    dist_from_defect = np.linalg.norm(_minimum_image(ref_pos - ref_pos[defect_atom], np.array([L, L, L])), axis=1)
    near = dist_from_defect < 4.0
    far = dist_from_defect > 12.0

    print(f"Single-atom defect kick (|kick| = {np.linalg.norm(kick):.2f} Ang):")
    print(f"  D2_min at/near the defect (< 4 Ang, n={near.sum()}): mean = {np.nanmean(d2_defect[near]):.3f}")
    print(f"  D2_min far from the defect (> 12 Ang, n={far.sum()}): mean = {np.nanmean(d2_defect[far]):.2e}")
    assert np.nanmean(d2_defect[near]) > 100 * np.nanmean(d2_defect[far]), \
        "FAIL: D2_min should be sharply localized near the defect"
    print("PASS: D2_min is sharply localized at the site of the actual non-affine rearrangement, "
          "near-zero elsewhere.\n")

    print("All self-tests passed.\n")

    # Test 3: same defect-kick test, but with a NONZERO box origin --
    # this specific case is what caught the real bug on first use
    # against run_1 (cKDTree periodic mode requires coordinates in
    # [0, boxsize), which raw LAMMPS coordinates essentially never
    # satisfy, since LAMMPS box origins are rarely exactly zero). The
    # first two tests above used an origin-anchored box and would NOT
    # have caught this, which is exactly why this third test exists.
    print("Self-test 3: same defect-kick check, with a nonzero box origin (as in real LAMMPS dumps).")
    origin = np.array([137.4, -22.9, 8.1])  # arbitrary, deliberately not zero and not simply related to L
    ref_pos_off = ref_pos + origin
    ref_box_off = ref_box + origin[:, None]
    cur_pos_off = cur_pos_defect + origin
    cur_box_off = cur_box_affine + origin[:, None]

    _, d2_off = compute_d2min_frame(ids, ref_pos_off, ref_box_off, ids, cur_pos_off, cur_box_off, cutoff=4.0)
    print(f"  D2_min at/near the defect (< 4 Ang): mean = {np.nanmean(d2_off[near]):.3f} "
          f"(should match Test 2's 0.939)")
    print(f"  D2_min far from the defect (> 12 Ang): mean = {np.nanmean(d2_off[far]):.2e} (should match Test 2's ~0)")
    assert np.allclose(d2_off[near], d2_defect[near], equal_nan=True), \
        "FAIL: result should be identical to Test 2 regardless of box origin"
    assert np.allclose(d2_off[far], d2_defect[far], equal_nan=True)
    print("PASS: identical result with a nonzero box origin -- this is the real-dump case.\n")

    print("All self-tests passed.")
