"""
find_seed_centers.py

Locates the inserted alpha-Mg crystal seed(s) in a composite's FIRST
tension-trajectory frame (strain=0, before any loading), using only
tools already built and validated in this project -- no OVITO
required. This unblocks interphase_shell_analysis.py, which needs a
seed center and radius as input.

Two independent, combined signals distinguish the crystal from the
amorphous matrix:
  1. q6 bond-orientational order close to the validated bulk HCP
     reference (~0.48-0.58; see bond_order_params.py) -- the matrix is
     amorphous and shows much lower, noisier q6.
  2. Local Mg purity: the seed is pure alpha-Mg (Methods 2.3), while
     the surrounding matrix is a mixed Mg-Zn-Ca alloy, so seed atoms
     should have an unusually Mg-rich local neighborhood compared to
     the bulk matrix composition.
Atoms satisfying both are clustered spatially with DBSCAN (the same
clustering approach already used for shear-band and D2min cluster
detection elsewhere in this project) to separate distinct seeds in
multi-seed (volume-fraction) composites, and each cluster's centroid
and radius are reported.

Seed placement enforces a minimum distance from the periodic boundary
(Methods 2.3), so no real seed should ever straddle a periodic wrap --
this implementation assumes that and does not attempt periodic-aware
clustering, which would otherwise be a real complication.

Usage:
    python3 find_seed_centers.py vf_1/tension_trajectory.dump
"""

import sys
import numpy as np
from sklearn.cluster import DBSCAN
from bond_order_params import compute_q4_q6
from multiframe_solute_shap import build_coordination, MG_TYPE, ZN_TYPE, CA_TYPE


def parse_first_frame(dump_path):
    with open(dump_path) as f:
        assert f.readline().startswith("ITEM: TIMESTEP")
        timestep = int(f.readline().strip())
        f.readline()
        n_atoms = int(f.readline().strip())
        f.readline()
        box = np.array([list(map(float, f.readline().split())) for _ in range(3)])
        header = f.readline().strip().replace("ITEM: ATOMS", "").split()
        idx = {c: header.index(c) for c in ("id", "type", "x", "y", "z")}
        types = np.empty(n_atoms, dtype=np.int64)
        pos = np.empty((n_atoms, 3))
        for i in range(n_atoms):
            parts = f.readline().split()
            types[i] = int(parts[idx["type"]])
            pos[i] = [float(parts[idx["x"]]), float(parts[idx["y"]]), float(parts[idx["z"]])]
    return timestep, box, types, pos


def find_seeds(dump_path, q6_threshold=0.40, mg_purity_threshold=0.95,
               cutoff=4.0, dbscan_eps=6.0, dbscan_min_samples=15):
    # eps/min_samples tuned via direct diagnostic on synthetic data: the
    # candidate-atom point cloud has real gaps near the seed surface
    # (median nearest-neighbor spacing ~3.2 Ang, but up to ~17 Ang at the
    # sparsest points), so eps=4.0 (tried first) failed to connect enough
    # of the cloud into one cluster; eps=6.0 with min_samples=15 gave
    # stable, consistent clustering across a range of parameter choices.
    timestep, box, types, pos = parse_first_frame(dump_path)

    q4, q6 = compute_q4_q6(pos, box, cutoff)
    n_mg, n_zn, n_ca = build_coordination(types, pos, box, cutoff)
    n_total = n_mg + n_zn + n_ca
    with np.errstate(invalid="ignore", divide="ignore"):
        mg_purity = n_mg / n_total

    is_mg_type = types == MG_TYPE
    candidate = is_mg_type & (q6 > q6_threshold) & (mg_purity > mg_purity_threshold) & (n_total > 0)
    print(f"Candidate crystal atoms: {candidate.sum()} of {len(types)} total "
          f"({is_mg_type.sum()} are Mg-type)")

    if candidate.sum() < dbscan_min_samples:
        print("No seed found (fewer candidate atoms than dbscan_min_samples) -- "
              "this may be a single-seed, near-zero-VF composition, or thresholds may need adjusting.")
        return []

    coords = pos[candidate]
    labels = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples).fit_predict(coords)
    seeds = []
    for lab in sorted(set(labels) - {-1}):
        cluster = coords[labels == lab]
        center = cluster.mean(axis=0)
        radius = np.percentile(np.linalg.norm(cluster - center, axis=1), 90)
        seeds.append(dict(center=tuple(center), radius=float(radius), n_atoms=len(cluster)))
        print(f"  Seed found: center=({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}), "
              f"radius~{radius:.2f} Ang, n_atoms={len(cluster)}")
    return seeds


if __name__ == "__main__":
    if len(sys.argv) > 1:
        find_seeds(sys.argv[1])
    else:
        print("Self-test: synthetic composite with a known HCP seed embedded in an amorphous-like matrix.\n")
        rng = np.random.default_rng(23)
        n_matrix = 8000
        L = 60.0
        origin = np.array([15.0, -20.0, 5.0])
        matrix_pos = rng.uniform(0, L, size=(n_matrix, 3)) + origin
        matrix_types = rng.choice([MG_TYPE, ZN_TYPE, CA_TYPE], size=n_matrix, p=[0.60, 0.35, 0.05])

        # Build a real (not fake) HCP lattice region, reusing the same
        # validated basis from the elastic-constants and bond-order work,
        # placed well clear of the box boundary per the real placement rule.
        a0, ca_ratio = 3.209, 1.619
        a1 = np.array([1.0, 0.0, 0.0]) * a0
        a2 = np.array([0.0, np.sqrt(3), 0.0]) * a0
        a3 = np.array([0.0, 0.0, ca_ratio]) * a0
        hcp_basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 5/6, 0.5), (0, 1/3, 0.5)]
        true_center = origin + np.array([30.0, 30.0, 30.0])
        true_radius = 12.0
        seed_pos = []
        for i in range(-5, 6):
            for j in range(-5, 6):
                for k in range(-5, 6):
                    base = i * a1 + j * a2 + k * a3
                    for f in hcp_basis:
                        p = base + f[0] * a1 + f[1] * a2 + f[2] * a3
                        if np.linalg.norm(p) < true_radius:
                            seed_pos.append(p + true_center)
        seed_pos = np.array(seed_pos)
        seed_types = np.full(len(seed_pos), MG_TYPE)
        print(f"Synthetic seed: {len(seed_pos)} atoms, true center={tuple(true_center.round(2))}, "
              f"true radius={true_radius}")

        # Remove matrix atoms that would overlap the inserted seed (matching
        # the real construction process's own overlap-removal step).
        from scipy.spatial import cKDTree
        tree = cKDTree(seed_pos)
        d, _ = tree.query(matrix_pos, k=1)
        keep = d > 2.0
        all_pos = np.vstack([matrix_pos[keep], seed_pos])
        all_types = np.concatenate([matrix_types[keep], seed_types])

        box = np.column_stack([origin, origin + L])
        tmp = "synthetic_composite.dump"
        with open(tmp, "w") as f:
            f.write("ITEM: TIMESTEP\n0\n")
            f.write("ITEM: NUMBER OF ATOMS\n"); f.write(f"{len(all_pos)}\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"{box[0,0]} {box[0,1]}\n{box[1,0]} {box[1,1]}\n{box[2,0]} {box[2,1]}\n")
            f.write("ITEM: ATOMS id type x y z\n")
            for i in range(len(all_pos)):
                f.write(f"{i+1} {all_types[i]} {all_pos[i,0]:.5f} {all_pos[i,1]:.5f} {all_pos[i,2]:.5f}\n")

        found = find_seeds(tmp)
        import os
        os.remove(tmp)

        assert len(found) == 1, f"FAIL: expected exactly 1 seed found, got {len(found)}"
        found_center = np.array(found[0]["center"])
        center_error = np.linalg.norm(found_center - true_center)
        radius_error = abs(found[0]["radius"] - true_radius)
        print(f"\nCenter error: {center_error:.2f} Ang, radius error: {radius_error:.2f} Ang")
        assert center_error < 2.0, f"FAIL: center error too large ({center_error:.2f} Ang)"
        assert radius_error < 3.0, f"FAIL: radius error too large ({radius_error:.2f} Ang)"
        print("PASS: seed correctly identified, with center and radius both close to ground truth.")
