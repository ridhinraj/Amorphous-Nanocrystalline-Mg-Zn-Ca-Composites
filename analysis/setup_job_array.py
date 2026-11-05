import os
import numpy as np

# ============================================================
# setup_job_array.py
#
# Generates two sets of run directories:
#
#   run_1 .. run_10   -- composition sweep, single small seed each
#                         (~0.7% VF), full 4-stage pipeline
#                         (quench -> build -> relax -> tension)
#
#   vf_1 .. vf_3      -- volume-fraction sweep on the BASELINE
#                         composition only (run_1's composition),
#                         reusing run_1's finished quenched glass
#                         instead of re-quenching. Only 3 stages
#                         (build -> relax -> tension) -- Stage 1
#                         is ~95% of total runtime, so skipping it
#                         for these 3 runs is both cheaper and
#                         scientifically cleaner (isolates VF as
#                         the only varying factor).
#
# IMPORTANT: vf_1..vf_3 depend on run_1/quenched_matrix_300K.data
# already existing. Do not submit submit_vf_array.sh until run_1's
# melt-quench stage has finished. See submit_vf_array.sh.
# ============================================================

# ---------------------------------------------------------
# PART A -- composition matrix (unchanged from earlier version)
# ---------------------------------------------------------
# Format: [xZn, xCa, Seed_Radius]
composition_parameters = [
    [0.35, 0.05, 25.0], # Run 1 (Baseline -- also the source for vf_1..vf_3)
    [0.30, 0.04, 20.0], # Run 2
    [0.40, 0.06, 30.0], # Run 3
    [0.25, 0.05, 25.0], # Run 4
    [0.35, 0.03, 15.0], # Run 5
    [0.35, 0.07, 35.0], # Run 6
    [0.30, 0.06, 25.0], # Run 7
    [0.40, 0.04, 20.0], # Run 8
    [0.32, 0.05, 28.0], # Run 9
    [0.38, 0.05, 22.0], # Run 10
]

# ---------------------------------------------------------
# PART B -- volume-fraction targets for the baseline composition
# ---------------------------------------------------------
# Deliberately modest (5/10/15%), not a full 0-40% sweep:
#  - continuum percolation threshold for randomly placed spheres
#    is around 29% volume fraction -- above that, seeds start
#    networking into a bicontinuous structure rather than staying
#    dispersed, which would change the deformation mechanism the
#    study is about.
#  - the packing algorithm below is plain rejection sampling,
#    which gets slow and unreliable as it approaches typical
#    random-sequential-packing saturation (~35-40% of the
#    available placement volume); 15% keeps comfortably clear
#    of that.
vf_targets = [0.05, 0.10, 0.15]
vf_seed_radius = 12.0   # smaller, fixed seed size so enough seeds fit to hit these VFs

# ============================================================
# TEMPLATES (token-substituted, not f-strings -- the builder
# script below has enough native Python braces/dict-comprehensions
# that f-string escaping would be error-prone; plain .replace()
# on unique tokens is safer).
# ============================================================

MQ_TEMPLATE = """include         ../common.mod
variable        Tmelt    equal 1500.0
variable        Tfinal   equal 300.0
variable        coolrate equal 1.0e12
variable        xZn      equal __XZN__
variable        xCa      equal __XCA__

lattice         fcc 4.5
region          box block 0 46 0 46 0 46 units lattice
create_box      3 box

include         ../potential.mod

create_atoms    1 box
set             type 1 type/fraction 2 ${xZn} 12345
set             type 1 type/fraction 3 $(v_xCa/(1.0-v_xZn)) 54321

minimize        1.0e-10 1.0e-10 10000 10000

# Frequent alternating pair: crash-recovery safety net only.
restart         50000 mq_restart.a mq_restart.b

velocity        all create 300.0 9834 dist gaussian
fix             1 all npt temp 300.0 ${Tmelt} 0.1 iso 0.0 0.0 5.0
timestep        0.0005
thermo          1000
thermo_style    custom step temp pe press vol density
run             100000
unfix           1
write_restart   ckpt_1_heated.restart

fix             2 all npt temp ${Tmelt} ${Tmelt} 0.1 iso 0.0 0.0 5.0
run             200000
unfix           2
write_restart   ckpt_2_melt_equil.restart

# NOTE: this fix is NPT, not NVT. NVT here was the bug found
# mid-project -- locking box volume during cooling forces the
# contracting liquid into an artificially tensile, under-density
# state as temperature drops. NPT lets volume track density.
variable        dtcool_ps   equal "(v_Tmelt - v_Tfinal) / v_coolrate * 1.0e12"
variable        nsteps_cool equal "round(v_dtcool_ps / 0.0005)"
fix             3 all npt temp ${Tmelt} ${Tfinal} 0.1 iso 0.0 0.0 5.0
run             ${nsteps_cool}
unfix           3
write_restart   ckpt_3_quenched.restart

fix             4 all npt temp ${Tfinal} ${Tfinal} 0.1 iso 0.0 0.0 5.0
run             100000
unfix           4
write_restart   ckpt_4_300K_settle.restart

fix             5 all npt temp 300.0 600.0 0.1 iso 0.0 0.0 5.0
run             200000
unfix           5

fix             6 all npt temp 600.0 600.0 0.1 iso 0.0 0.0 5.0
run             1000000
unfix           6
write_restart   ckpt_6_600K_anneal.restart

fix             7 all npt temp 600.0 300.0 0.1 iso 0.0 0.0 5.0
run             600000
unfix           7

fix             8 all npt temp 300.0 300.0 0.1 iso 0.0 0.0 5.0
run             200000
unfix           8

write_data      quenched_matrix_300K.data
"""

BUILDER_TEMPLATE = '''import numpy as np
from scipy.spatial import cKDTree
from pymatgen.core import Structure, Lattice
from pymatgen.io.lammps.data import LammpsData
import pandas as pd

AMORPHOUS_FILE = "__AMORPHOUS_FILE__"
SEED_RADIUS    = __SEED_RADIUS__
TARGET_VF      = __TARGET_VF__   # None -> single seed at box center; float -> packed multi-seed
MIN_GAP        = 4.0
WALL_MARGIN    = 4.0
OVERLAP_CUTOFF = 2.0
RNG_SEED       = __RNG_SEED__


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
'''

RELAX_TEMPLATE = """include         ../common.mod
read_data       MgZnCa_Composite_Ready.data
include         ../potential.mod

reset_timestep  0
thermo          100
thermo_style    custom step temp pe press vol density
minimize        1.0e-10 1.0e-10 10000 10000

timestep        0.001
velocity        all create 300.0 12345 rot yes dist gaussian
fix             1 all npt temp 300.0 300.0 0.1 iso 0.0 0.0 5.0
run             50000

write_data      AN_composite_relaxed_300K.data
"""

TENSION_TEMPLATE = """include         ../common.mod
read_data       AN_composite_relaxed_300K.data
include         ../potential.mod

variable        temp equal 300.0
variable        srate equal 1.0e9
variable        erate equal v_srate/1.0e12

reset_timestep  0
timestep        0.001
velocity        all create ${temp} 87287 dist gaussian

fix             1 all npt temp ${temp} ${temp} 0.1 y 0.0 0.0 5.0 z 0.0 0.0 5.0
run             20000
unfix           1

# fix 2 keeps y/z pressure-coupled through the pulling run (uniaxial
# stress, free lateral contraction) instead of switching to plain nvt
# (which would freeze y/z and impose uniaxial strain instead).
fix             2 all npt temp ${temp} ${temp} 0.1 y 0.0 0.0 5.0 z 0.0 0.0 5.0
fix             3 all deform 1 x erate ${erate} units box remap x

# L0 is frozen at the box length *now* via immediate evaluation
# $(lx) -- "variable L0 equal lx" (without $()) would silently
# re-evaluate to the current lx every time it's referenced, making
# strain identically zero for the whole run.
variable        L0 equal $(lx)
variable        strain equal "(lx - v_L0)/v_L0"
variable        stress_x equal "-pxx/10000"
variable        stress_y equal "-pyy/10000"
variable        stress_z equal "-pzz/10000"

compute         peratom_stress all stress/atom NULL
thermo          1000
thermo_style    custom step v_strain v_stress_x v_stress_y v_stress_z temp pe press

fix             4 all print 100 "${strain} ${stress_x} ${stress_y} ${stress_z}" file tension_stress_strain.txt screen no title "Strain Stress_X_GPa Stress_Y_GPa Stress_Z_GPa"
dump            1 all custom 5000 tension_trajectory.dump id type x y z c_peratom_stress[1] c_peratom_stress[2] c_peratom_stress[3] c_peratom_stress[4] c_peratom_stress[5] c_peratom_stress[6]

run             150000
"""

# ============================================================
# GENERATE PART A: run_1 .. run_10
# ============================================================
for i, (xZn, xCa, radius) in enumerate(composition_parameters):
    run_id = i + 1
    run_dir = f"run_{run_id}"
    os.makedirs(run_dir, exist_ok=True)

    mq_script = MQ_TEMPLATE.replace("__XZN__", str(xZn)).replace("__XCA__", str(xCa))
    with open(f"{run_dir}/week2_melt_quench.in", "w") as f:
        f.write(mq_script)

    builder_script = (
        BUILDER_TEMPLATE
        .replace("__AMORPHOUS_FILE__", "quenched_matrix_300K.data")
        .replace("__SEED_RADIUS__", str(radius))
        .replace("__TARGET_VF__", "None")
        .replace("__RNG_SEED__", str(1000 + run_id))
    )
    with open(f"{run_dir}/build_an_composite.py", "w") as f:
        f.write(builder_script)

    with open(f"{run_dir}/in.relax_composite", "w") as f:
        f.write(RELAX_TEMPLATE)

    with open(f"{run_dir}/week4_tension.in", "w") as f:
        f.write(TENSION_TEMPLATE)

# ============================================================
# GENERATE PART B: vf_1 .. vf_3 (baseline composition, no re-quench)
# ============================================================
for i, target_vf in enumerate(vf_targets):
    run_id = i + 1
    run_dir = f"vf_{run_id}"
    os.makedirs(run_dir, exist_ok=True)

    builder_script = (
        BUILDER_TEMPLATE
        .replace("__AMORPHOUS_FILE__", "../run_1/quenched_matrix_300K.data")
        .replace("__SEED_RADIUS__", str(vf_seed_radius))
        .replace("__TARGET_VF__", str(target_vf))
        .replace("__RNG_SEED__", str(2000 + run_id))
    )
    with open(f"{run_dir}/build_an_composite.py", "w") as f:
        f.write(builder_script)

    with open(f"{run_dir}/in.relax_composite", "w") as f:
        f.write(RELAX_TEMPLATE)

    with open(f"{run_dir}/week4_tension.in", "w") as f:
        f.write(TENSION_TEMPLATE)

print("Generated run_1..run_10 (composition sweep) and vf_1..vf_3 (volume-fraction sweep).")
print("vf_1..vf_3 require run_1/quenched_matrix_300K.data to exist before they can run.")
