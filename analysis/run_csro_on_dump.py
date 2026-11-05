"""
run_csro_on_dump.py

Computes real global (bulk) and per-atom Warren-Cowley CSRO
parameters at a given strain, on an actual tension_trajectory.dump.
Reuses the exact same frame-parsing function already used by
run_bond_order_on_dump.py, so both tools read dumps identically.

Usage:
    python3 run_csro_on_dump.py run_1/tension_trajectory.dump 0.071 --erate_per_ps 0.001 --cutoff 4.0
"""
import argparse
import numpy as np
from run_bond_order_on_dump import parse_dump_frame_at_strain
from warren_cowley_csro import global_wc_parameters, per_atom_own_type_wc, TYPE_NAMES
from multiframe_solute_shap import ZN_TYPE, CA_TYPE

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dump_path")
    parser.add_argument("target_strain", type=float)
    parser.add_argument("--erate_per_ps", type=float, default=0.001)
    parser.add_argument("--cutoff", type=float, default=4.0,
                         help="neighbor-shell cutoff; 4.0 Ang comfortably covers the Ca-Ca RDF "
                              "peak (3.825 Ang) as well as Zn-Zn/Mg-Mg, avoiding the selection "
                              "bias a 3.25 Ang cutoff caused for Ca specifically (found directly "
                              "on run_1: a 3.25 Ang cutoff dropped 37% of Ca atoms for having "
                              "too few neighbors, versus <1% of Zn atoms).")
    args = parser.parse_args()

    (timestep, box, types, pos), actual_strain = parse_dump_frame_at_strain(
        args.dump_path, args.target_strain, args.erate_per_ps)
    print(f"Target strain {args.target_strain:.3f} -> nearest frame at strain {actual_strain:.3f} "
          f"(timestep {timestep})\n")

    print("=== Global (bulk) Warren-Cowley parameters ===")
    alpha, x = global_wc_parameters(types, pos, box, cutoff=args.cutoff)
    print(f"Bulk mole fractions: Mg={x[1]:.4f}, Zn={x[2]:.4f}, Ca={x[3]:.4f}\n")
    for pair, val in alpha.items():
        # alpha = 1 - P_AB/x_B: NEGATIVE means A-B pairs are MORE common
        # than random (A and B prefer each other -- ordering/intermixing);
        # POSITIVE means A-B pairs are LESS common than random (A and B
        # avoid each other). For a SAME-type pair (A=B) specifically,
        # positive therefore means that species avoids ITS OWN kind
        # (dispersion), and negative means it prefers its own kind
        # (clustering/segregation) -- the earlier version of this label
        # had this backwards for the same-type case.
        flag = ""
        if abs(val) > 0.05:
            flag = "  (favored: A-B pairs more common than random)" if val < 0 else \
                   "  (avoided: A-B pairs less common than random)"
        print(f"  alpha_{pair[0]}-{pair[1]:<3} = {val:+.4f}{flag}")

    print("\n=== Per-atom own-type Warren-Cowley parameter (Ca vs Zn) ===")
    per_atom_alpha = per_atom_own_type_wc(types, pos, box, cutoff=args.cutoff)
    is_zn = (types == ZN_TYPE) & ~np.isnan(per_atom_alpha)
    is_ca = (types == CA_TYPE) & ~np.isnan(per_atom_alpha)
    print(f"n_Zn={is_zn.sum()}, n_Ca={is_ca.sum()}")
    print(f"  mean alpha_Zn-Zn = {np.mean(per_atom_alpha[is_zn]):+.4f}  "
          f"(std {np.std(per_atom_alpha[is_zn]):.4f})")
    print(f"  mean alpha_Ca-Ca = {np.mean(per_atom_alpha[is_ca]):+.4f}  "
          f"(std {np.std(per_atom_alpha[is_ca]):.4f})")

    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(per_atom_alpha[is_ca], per_atom_alpha[is_zn])
    print(f"  Ca - Zn per-atom alpha difference: "
          f"{np.mean(per_atom_alpha[is_ca]) - np.mean(per_atom_alpha[is_zn]):+.4f}  "
          f"(Mann-Whitney p={p:.2e})")
