"""
compute_rdf.py

Computes the total (and optionally partial) radial distribution
function g(r) for a quenched LAMMPS structure, with periodic boundary
conditions handled correctly via scipy's cKDTree boxsize support
(confirmed working directly before using it here -- boxsize periodic
wrapping is not obvious from cKDTree's signature via introspection,
so it was tested against a small synthetic case first).

This is the piece Section 3.2 (amorphous structure characterization)
actually needs and doesn't have yet -- nothing in this project has
computed a real g(r) curve before now. Run this against each
composition's quenched_matrix_300K.data (the Stage 1 output, before
composite building) to get the real data Figure 2 needs.

Usage:
    python3 compute_rdf.py run_1/quenched_matrix_300K.data
    python3 compute_rdf.py run_1/quenched_matrix_300K.data --partial
"""

import sys
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

TYPE_NAMES = {1: "Mg", 2: "Zn", 3: "Ca"}


def read_lammps_data(path):
    """
    Minimal, dependency-free parser for a LAMMPS 'data' file (the
    format written by write_data -- quenched_matrix_300K.data,
    AN_composite_relaxed_300K.data, etc.). Deliberately avoids
    pymatgen here: this project already has enough environment
    friction between the anaconda module (has pymatgen) and shap_env
    (doesn't, and adding it risks reopening the numpy-version
    conflict between the two). Only reads what compute_gr actually
    needs -- box bounds and atom positions/types -- and skips
    Masses/Velocities/etc. sections entirely.
    """
    with open(path) as f:
        lines = f.readlines()

    n_atoms = None
    box_bounds = {}
    atoms_start = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.endswith("atoms") and n_atoms is None:
            n_atoms = int(stripped.split()[0])
        elif "xlo xhi" in stripped:
            parts = stripped.split()
            box_bounds["x"] = (float(parts[0]), float(parts[1]))
        elif "ylo yhi" in stripped:
            parts = stripped.split()
            box_bounds["y"] = (float(parts[0]), float(parts[1]))
        elif "zlo zhi" in stripped:
            parts = stripped.split()
            box_bounds["z"] = (float(parts[0]), float(parts[1]))
        elif stripped.startswith("Atoms"):
            atoms_start = i + 2  # skip the "Atoms" line and the blank line after it
            break

    if n_atoms is None or atoms_start is None or not box_bounds:
        raise ValueError(f"Could not parse header of {path} -- missing atom count, "
                          f"box bounds, or Atoms section marker.")

    positions = np.zeros((n_atoms, 3))
    types = np.zeros(n_atoms, dtype=int)
    for k in range(n_atoms):
        parts = lines[atoms_start + k].split()
        # Standard atomic-style: id type x y z [extra columns ignored]
        types[k] = int(parts[1])
        positions[k] = [float(parts[2]), float(parts[3]), float(parts[4])]

    box_length_x = box_bounds["x"][1] - box_bounds["x"][0]
    box_length_y = box_bounds["y"][1] - box_bounds["y"][0]
    box_length_z = box_bounds["z"][1] - box_bounds["z"][0]

    # Shift so the box origin is at (0,0,0) -- scipy's cKDTree periodic
    # boxsize wrapping requires coordinates in [0, boxsize), but real
    # LAMMPS boxes almost never start at exactly 0 (e.g. xlo=1.64 in
    # this project's actual output, an NPT box-relaxation artifact,
    # not an error). Missed in the original synthetic tests because
    # every one of them happened to use a box starting at 0.0 --
    # caught only once run against real data with a real, nonzero
    # origin. Modulo (not a plain subtraction) handles any atom that
    # ends up sitting exactly on the upper boundary after shifting.
    origin = np.array([box_bounds["x"][0], box_bounds["y"][0], box_bounds["z"][0]])
    box_lengths_arr = np.array([box_length_x, box_length_y, box_length_z])
    positions = (positions - origin) % box_lengths_arr

    return positions, types, (box_length_x, box_length_y, box_length_z)


def compute_gr(positions, box_length, r_max=10.0, n_bins=200, type_pair=None, types=None):
    """
    Standard g(r) normalization: g(r) = 2*n_pairs(r,r+dr) / (N * rho * 4*pi*r^2*dr),
    where n_pairs is the count of UNIQUE pairs in that shell (query_pairs
    returns each pair once, hence the factor of 2 -- each pair
    contributes to both particles' local coordination shells).

    If type_pair=(a, b) is given, restricts to distances between atoms
    of type a and type b specifically (partial RDF); otherwise computes
    the total RDF over all atoms.
    """
    N = len(positions)
    V = box_length ** 3
    rho = N / V

    tree = cKDTree(positions, boxsize=box_length)
    pairs = tree.query_pairs(r=r_max, output_type="ndarray")

    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    if type_pair is not None:
        a, b = type_pair
        mask = ((types[i_idx] == a) & (types[j_idx] == b)) | ((types[i_idx] == b) & (types[j_idx] == a))
        i_idx, j_idx = i_idx[mask], j_idx[mask]

    delta = positions[i_idx] - positions[j_idx]
    delta -= box_length * np.round(delta / box_length)  # minimum image
    dists = np.linalg.norm(delta, axis=1)

    bin_edges = np.linspace(0, r_max, n_bins + 1)
    counts, _ = np.histogram(dists, bins=bin_edges)

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dr = bin_edges[1] - bin_edges[0]
    shell_volume = 4 * np.pi * bin_centers ** 2 * dr

    if type_pair is not None:
        # partial RDF: normalize by the density of the SECOND type,
        # per standard partial-RDF convention, using N of the relevant
        # type pair for the "per reference atom" normalization
        n_a = np.sum(types == type_pair[0])
        n_b = np.sum(types == type_pair[1])
        if type_pair[0] == type_pair[1]:
            g_r = 2 * counts / (n_a * (n_b / V) * shell_volume)
        else:
            g_r = counts / (n_a * (n_b / V) * shell_volume)
    else:
        g_r = 2 * counts / (N * rho * shell_volume)

    return bin_centers, g_r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file")
    parser.add_argument("--r_max", type=float, default=10.0)
    parser.add_argument("--n_bins", type=int, default=200)
    parser.add_argument("--partial", action="store_true", help="also compute Mg-Mg, Zn-Zn, Ca-Ca, Mg-Zn, Mg-Ca, Zn-Ca partials")
    parser.add_argument("--out", default=None, help="output CSV path (default: <input>_rdf.csv)")
    args = parser.parse_args()

    print(f"Loading {args.data_file} ...")
    positions, types, box_lengths = read_lammps_data(args.data_file)
    box_length = box_lengths[0]

    if not (np.isclose(box_lengths[0], box_lengths[1], rtol=1e-6) and
            np.isclose(box_lengths[0], box_lengths[2], rtol=1e-6)):
        print(f"WARNING: box is not cubic to high precision (x={box_lengths[0]:.4f}, "
              f"y={box_lengths[1]:.4f}, z={box_lengths[2]:.4f}) -- using x-dimension length.")

    print(f"N = {len(positions)}, box length = {box_length:.4f} Ang, "
          f"density = {len(positions) / box_length**3:.6f} atoms/Ang^3")

    r, g_total = compute_gr(positions, box_length, args.r_max, args.n_bins)

    results = {"r": r, "g_total": g_total}

    if args.partial:
        pairs_to_compute = [(1, 1, "Mg-Mg"), (2, 2, "Zn-Zn"), (3, 3, "Ca-Ca"),
                             (1, 2, "Mg-Zn"), (1, 3, "Mg-Ca"), (2, 3, "Zn-Ca")]
        for a, b, label in pairs_to_compute:
            if np.sum(types == a) > 0 and np.sum(types == b) > 0:
                _, g_partial = compute_gr(positions, box_length, args.r_max, args.n_bins,
                                           type_pair=(a, b), types=types)
                results[f"g_{label}"] = g_partial

    df = pd.DataFrame(results)
    out_path = args.out or args.data_file.replace(".data", "_rdf.csv")
    df.to_csv(out_path, index=False)

    first_peak_idx = np.argmax(g_total[r > 1.5])  # skip the unphysical r->0 region
    first_peak_r = r[r > 1.5][first_peak_idx]
    first_peak_g = g_total[r > 1.5][first_peak_idx]
    print(f"\nFirst peak: r = {first_peak_r:.3f} Ang, g(r) = {first_peak_g:.3f}")
    print(f"g(r) at r_max ({args.r_max} Ang) = {g_total[-1]:.3f} (should be close to 1.0 for a converged, "
          f"correctly-normalized RDF at this range)")
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
