import numpy as np
from scipy.spatial import cKDTree
from pymatgen.core import Structure, Lattice
from pymatgen.io.lammps.data import LammpsData
import pandas as pd

AMORPHOUS_FILE = "../run_5/quenched_matrix_300K.data"
SEED_RADIUS    = 12.0
TARGET_VF      = 0.15   # None -> single seed at box center; float -> packed multi-seed
MIN_GAP        = 4.0
WALL_MARGIN    = 4.0
OVERLAP_CUTOFF = 2.0
RNG_SEED       = 4003


def create_hcp_mg_seed(radius):
    a, c = 3.20, 5.21
    lattice = Lattice.hexagonal(a, c)
    mg_unit = Structure(lattice, ["Mg", "Mg"], [[0, 0, 0], [1 / 3, 2 / 3, 1 / 2]])
    n_a = int(np.ceil(2 * radius / a)) + 2
    n_c = int(np.ceil(2 * radius / c)) + 2
    mg_supercell = mg_unit.make_supercell([n_a, n_a, n_c], in_place=False)
    coords = mg_supercell.cart_coords
    center = np.mean(coords, axis=0)
    shifted_coords = coords - center
    distances = np.linalg.norm(shifted_coords, axis=1)
    return shifted_coords[np.where(distances <= radius)[0]]


def generate_seed_centers(box_bounds, target_vf, seed_radius, min_gap, wall_margin, rng):
    """
    Rejection-sample non-overlapping sphere centers to approximately
    hit target_vf inside the actual simulation box (read from the
    LAMMPS data file, not assumed). Centers are kept at least
    (seed_radius + wall_margin) from every box face so seeds never
    need to wrap across the periodic boundary -- trades a little
    usable volume for not having to split a seed's atoms across a
    periodic image, which is fine at the modest VF targets used here.
    """
    box_bounds = np.asarray(box_bounds)  # shape (3,2): [[xlo,xhi],[ylo,yhi],[zlo,zhi]]
    box_volume = np.prod(box_bounds[:, 1] - box_bounds[:, 0])
    seed_volume = (4.0 / 3.0) * np.pi * seed_radius ** 3
    n_target = max(1, int(round(target_vf * box_volume / seed_volume)))

    lo = box_bounds[:, 0] + seed_radius + wall_margin
    hi = box_bounds[:, 1] - seed_radius - wall_margin
    if np.any(hi <= lo):
        raise ValueError("seed_radius too large for this box at this wall_margin")

    centers = []
    min_center_dist = 2 * seed_radius + min_gap
    max_attempts = 200000
    attempts = 0
    while len(centers) < n_target and attempts < max_attempts:
        attempts += 1
        candidate = rng.uniform(lo, hi)
        ok = True
        for c in centers:
            if np.linalg.norm(candidate - c) < min_center_dist:
                ok = False
                break
        if ok:
            centers.append(candidate)

    achieved_vf = len(centers) * seed_volume / box_volume
    print(f"Requested VF={target_vf:.3f}, target seeds={n_target}, "
          f"placed={len(centers)} in {attempts} attempts, achieved VF={achieved_vf:.4f}")
    if len(centers) < n_target:
        print("WARNING: rejection sampling did not reach the target seed count "
              "before max_attempts -- achieved VF is below target. Reduce "
              "target_vf, reduce seed_radius, or raise max_attempts.")
    return np.array(centers)


def assemble_an_composite(amorph_data, seed_centers, seed_radius, overlap_cutoff):
    atoms_df = amorph_data.atoms
    amorph_coords = atoms_df[["x", "y", "z"]].values
    amorph_types = atoms_df["type"].values

    keep_mask = np.ones(len(amorph_coords), dtype=bool)
    for center in seed_centers:
        d = np.linalg.norm(amorph_coords - center, axis=1)
        keep_mask &= d > seed_radius
    matrix_coords = amorph_coords[keep_mask]
    matrix_types = amorph_types[keep_mask]

    all_seed_coords = [create_hcp_mg_seed(seed_radius) + center for center in seed_centers]
    seed_abs_coords = np.vstack(all_seed_coords)
    # Mg is mass type 1 in this mapping.
    seed_types = np.full(len(seed_abs_coords), 1, dtype=int)

    # Drop matrix atoms too close to any seed atom.
    tree_seed = cKDTree(seed_abs_coords)
    tree_matrix = cKDTree(matrix_coords)
    overlap_pairs = tree_matrix.query_ball_tree(tree_seed, r=overlap_cutoff)
    overlap_idx = [i for i, n in enumerate(overlap_pairs) if len(n) > 0]
    clean_mask = np.ones(len(matrix_coords), dtype=bool)
    clean_mask[overlap_idx] = False

    # Drop the rare seed-seed atom overlap at adjoining seed surfaces.
    # Center-distance rejection sampling prevents most of this, but
    # doesn't guarantee zero overlap at the atomic scale.
    seed_pairs = tree_seed.query_pairs(r=overlap_cutoff)
    seed_drop = {j for i, j in seed_pairs}
    seed_keep_mask = np.ones(len(seed_abs_coords), dtype=bool)
    for idx in seed_drop:
        seed_keep_mask[idx] = False

    final_coords = np.vstack([matrix_coords[clean_mask], seed_abs_coords[seed_keep_mask]])
    final_types = np.hstack([matrix_types[clean_mask], seed_types[seed_keep_mask]])

    # NOTE: pymatgen's LammpsData expects the atoms DataFrame to be
    # INDEXED by atom id (1, 2, 3, ...), with columns type/x/y/z --
    # not a separate "id" column with a default 0-based RangeIndex.
    # The earlier version here did the latter, which caused pymatgen
    # to write BOTH the stray default index and the "id" column as
    # separate fields, producing a 6-column Atoms section instead of
    # the 5 columns atom_style atomic requires. LAMMPS's read_data
    # then fails immediately with "Incorrect format in Atoms section"
    # the moment Stage 3 tries to read this file. Confirmed by
    # reproducing the exact bug and verifying this fix against a
    # real pymatgen round-trip before this was corrected.
    new_atoms_df = pd.DataFrame({
        "type": final_types,
        "x": final_coords[:, 0], "y": final_coords[:, 1], "z": final_coords[:, 2]
    }, index=pd.RangeIndex(start=1, stop=len(final_coords) + 1))

    amorph_data.atoms = new_atoms_df
    # The loaded quench data carries real MD velocities for its
    # original 389,344 atoms. After carving the cavity and inserting
    # the seed, the atom set (count and identity) is completely
    # different, so those velocities are stale and mismatched with
    # the new atom count -- write_file() would emit an inconsistent
    # Velocities section (wrong atom count vs. the header), which
    # LAMMPS/pymatgen both reject. They're harmless to drop: Stage 3
    # (in.relax_composite) re-initializes velocities from scratch
    # immediately after reading this file regardless.
    amorph_data.velocities = None
    amorph_data.write_file("MgZnCa_Composite_Ready.data")
    print(f"Wrote MgZnCa_Composite_Ready.data with {len(seed_centers)} seed(s), "
          f"{len(final_coords)} total atoms.")


if __name__ == "__main__":
    rng = np.random.default_rng(RNG_SEED)
    amorph_data = LammpsData.from_file(AMORPHOUS_FILE, atom_style="atomic")
    box_bounds = np.array(amorph_data.box.bounds)  # actual box, not assumed

    if TARGET_VF is None:
        box_center = box_bounds.mean(axis=1)
        centers = np.array([box_center])
    else:
        centers = generate_seed_centers(box_bounds, TARGET_VF, SEED_RADIUS, MIN_GAP, WALL_MARGIN, rng)

    assemble_an_composite(amorph_data, centers, SEED_RADIUS, OVERLAP_CUTOFF)
