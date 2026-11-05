"""
run_bond_order_on_dump.py

Computes real q4/q6 for every Zn and Ca atom at a given strain, on an
actual tension_trajectory.dump, and reports the Ca-vs-Zn comparison --
the same question as the existing SHAP analysis (does Ca sit in a
geometrically distinct local environment from Zn?), but from local
packing order instead of simple neighbor counts.

Usage:
    python3 run_bond_order_on_dump.py run_1/tension_trajectory.dump 0.071 --erate_per_ps 0.001
"""
import sys
import argparse
import numpy as np
from bond_order_params import compute_q4_q6

MG_TYPE, ZN_TYPE, CA_TYPE = 1, 2, 3
T0_STEPS = 20000
TIMESTEP_PS = 0.001


def frame_strain(timestep, erate_per_ps):
    return (timestep - T0_STEPS) * erate_per_ps * TIMESTEP_PS


def parse_dump_frame_at_strain(dump_path, target_strain, erate_per_ps):
    """Reads id, type, x, y, z for every frame, returns the one nearest
    target_strain (position-only -- no stress columns needed for this
    analysis)."""
    frames = []
    with open(dump_path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
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
            frames.append((timestep, box, types, pos))
    strains = np.array([frame_strain(fr[0], erate_per_ps) for fr in frames])
    i = int(np.argmin(np.abs(strains - target_strain)))
    return frames[i], strains[i]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dump_path")
    parser.add_argument("target_strain", type=float)
    parser.add_argument("--erate_per_ps", type=float, default=0.001)
    parser.add_argument("--cutoff", type=float, default=3.25,
                         help="first-neighbor-shell cutoff; 3.25 Ang matches this project's RDF-derived value")
    args = parser.parse_args()

    (timestep, box, types, pos), actual_strain = parse_dump_frame_at_strain(
        args.dump_path, args.target_strain, args.erate_per_ps)
    print(f"Target strain {args.target_strain:.3f} -> nearest frame at strain {actual_strain:.3f} "
          f"(timestep {timestep})")

    q4, q6 = compute_q4_q6(pos, box, args.cutoff)

    is_zn = (types == ZN_TYPE) & ~np.isnan(q4)
    is_ca = (types == CA_TYPE) & ~np.isnan(q4)
    print(f"\nn_Zn={is_zn.sum()}, n_Ca={is_ca.sum()}")
    print(f"{'':12}{'mean q4':>10}{'mean q6':>10}{'std q4':>10}{'std q6':>10}")
    print(f"{'Zn':<12}{np.nanmean(q4[is_zn]):>10.4f}{np.nanmean(q6[is_zn]):>10.4f}"
          f"{np.nanstd(q4[is_zn]):>10.4f}{np.nanstd(q6[is_zn]):>10.4f}")
    print(f"{'Ca':<12}{np.nanmean(q4[is_ca]):>10.4f}{np.nanmean(q6[is_ca]):>10.4f}"
          f"{np.nanstd(q4[is_ca]):>10.4f}{np.nanstd(q6[is_ca]):>10.4f}")

    from scipy.stats import mannwhitneyu
    for name, vals in [("q4", q4), ("q6", q6)]:
        u, p = mannwhitneyu(vals[is_ca], vals[is_zn])
        diff = np.nanmean(vals[is_ca]) - np.nanmean(vals[is_zn])
        print(f"Ca - Zn {name} difference: {diff:+.4f}  (Mann-Whitney p={p:.2e})")
