#!/usr/bin/env python3
"""
setup_new_md.py -- generate case directories for the additional MD runs that address the limitations
listed in Section 4.6 (confounded Ca series, single trajectories, one-point temperature test, thermal
history of Compositions 1/9, crystal-phase reference, cell size).

Run from the PROJECT ROOT on cluster 1 (the directory that holds common.mod, potential.mod, potentials/
and the run_* source directories).  It does NOT launch anything; it writes
    <case>/case.cfg  + the LAMMPS / builder inputs each case needs
    cases/<block>.txt         one case name per line (SLURM array lines)
    cases/<block>.ntasks      MPI ranks for that block
    cases_manifest.csv        what each case is for (used by collect_new_md.py)
Templates are the exact ones from setup_job_array.py (see md_templates.py), so every earlier fix carries over
(pd.RangeIndex, cleared velocities, frozen L0, y/z NPT coupling).
"""
import os, csv, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from md_templates import MQ_TEMPLATE, BUILDER_TEMPLATE, RELAX_TEMPLATE, TENSION_TEMPLATE

# composition table copied from setup_job_array.py: run -> (xZn, xCa, original single-seed radius)
COMP = {1: (0.35, 0.05, 25.0), 3: (0.40, 0.06, 30.0), 4: (0.25, 0.05, 25.0), 5: (0.35, 0.03, 15.0),
        6: (0.35, 0.07, 35.0), 9: (0.32, 0.05, 28.0)}
VF_R = 12.0            # matched seed radius used in the existing vf_1-3 and vf_r6_1-3 series

manifest, blocks = [], {}

def tension(temp=300.0, srate=1.0e9, equil=20000, strain_end=0.15):
    """The existing tension template with temperature / rate / equilibration / length changed. At the defaults it is
    identical to the one used for all 25 existing runs. 150,000 steps = 15 % strain at 1e9 /s with dt = 1 fs."""
    t = TENSION_TEMPLATE
    nsteps = round(strain_end / (srate * 1.0e-12 * 0.001))
    for a, b in [("variable        temp equal 300.0", f"variable        temp equal {temp}"),
                 ("variable        srate equal 1.0e9", f"variable        srate equal {srate:.1e}".replace("e+0", "e")),
                 ("run             20000", f"run             {equil}"),
                 ("run             150000", f"run             {nsteps}")]:
        assert t.count(a) == 1, a
        t = t.replace(a, b)
    return t

def builder(glass, radius, vf, rng):
    return (BUILDER_TEMPLATE.replace("__AMORPHOUS_FILE__", glass).replace("__SEED_RADIUS__", str(radius))
            .replace("__TARGET_VF__", "None" if vf is None else str(vf)).replace("__RNG_SEED__", str(rng)))

def quench(xzn, xca, seeds=None):
    t = MQ_TEMPLATE.replace("__XZN__", str(xzn)).replace("__XCA__", str(xca))
    if seeds:                                   # independent realization: new random assignment + new velocities
        s_zn, s_ca, s_v = seeds
        for a, b in [("12345", str(s_zn)), ("54321", str(s_ca)), ("9834", str(s_v))]:
            assert t.count(a) == 1; t = t.replace(a, b)
    return t

def write(case, name, text):
    os.makedirs(case, exist_ok=True)
    with open(f"{case}/{name}", "w") as f: f.write(text)

def add(block, case, stages, cfg, files, info, ntasks=40):
    cfg = dict(cfg, STAGES=" ".join(stages), BLOCK=block)
    write(case, "case.cfg", "".join(f'{k}="{v}"\n' for k, v in cfg.items()))
    for name, text in files.items(): write(case, name, text)
    blocks.setdefault(block, {"cases": [], "ntasks": ntasks})["cases"].append(case)
    manifest.append(dict(case=case, block=block, **info))

# --------------------------------------------------------------------------------------------------
# Block 1  ca_matched : xCa = 0.03 at the matched geometry (12 A seeds, 5/10/15 %), same as vf_1-3 / vf_r6_1-3.
#          Closes the missing 0.03 -> 0.05 matched-geometry control.   Source: run_5 glass (already quenched)
# --------------------------------------------------------------------------------------------------
for k, vf in enumerate([0.05, 0.10, 0.15], 1):
    c = f"vf_r5_{k}"
    add("ca_matched", c, ["build", "relax", "tension"],
        dict(SOURCE_GLASS="../run_5/quenched_matrix_300K.data", AXIS="x"),
        {"build_an_composite.py": builder("../run_5/quenched_matrix_300K.data", VF_R, vf, 4000 + k),
         "in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
        dict(group=f"r5_vf{int(vf*100)}", xZn=0.35, xCa=0.03, vf=vf, radius=VF_R, axis="x",
             purpose="matched-geometry Ca series (0.03)"))

# --------------------------------------------------------------------------------------------------
# Block 2  pack_reps : three further seed-placement realizations (new RNG seed) at 12 A / 10 % for each of the three Ca
#          levels, on the SAME quenched glass.  With the existing 10 % case that makes n = 4 per composition.
#          Measures placement scatter only (glass structure is shared) -- see README for what this does / does not show.
# --------------------------------------------------------------------------------------------------
for run in (5, 1, 6):
    for j, tag in enumerate("abc"):
        c = f"rep_r{run}_{tag}"
        add("pack_reps", c, ["build", "relax", "tension"],
            dict(SOURCE_GLASS=f"../run_{run}/quenched_matrix_300K.data", AXIS="x"),
            {"build_an_composite.py": builder(f"../run_{run}/quenched_matrix_300K.data", VF_R, 0.10, 5000 + run*10 + j),
             "in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
            dict(group=f"r{run}_vf10", xZn=COMP[run][0], xCa=COMP[run][1], vf=0.10, radius=VF_R, axis="x",
                 purpose="seed-placement replicate"))

# --------------------------------------------------------------------------------------------------
# Block 3  glass0 : 0 % volume fraction -- the quenched matrix itself, pulled with the identical protocol
# --------------------------------------------------------------------------------------------------
for run in (5, 1, 6):
    c = f"glass0_r{run}"
    add("glass0", c, ["copyglass", "relax", "tension"],
        dict(SOURCE_GLASS=f"../run_{run}/quenched_matrix_300K.data", AXIS="x"),
        {"in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
        dict(group=f"r{run}_vf0", xZn=COMP[run][0], xCa=COMP[run][1], vf=0.0, radius=0, axis="x",
             purpose="fully amorphous reference"))

# --------------------------------------------------------------------------------------------------
# Block 4  hcp : defect-free HCP Mg single-crystal cell (~339k atoms, ~20 nm), tension along [2-1-10] (x) and [0001] (z)
# --------------------------------------------------------------------------------------------------
HCP = """include         ../common.mod
# Defect-free HCP Mg, orthogonal 4-atom cell: lx = 62a, ly = 36*sqrt(3)a, lz = 38c  (~198 x 200 x 198 A, 339,264 atoms)
# Axis of tension: __AXIS__   (x = [2-1-10], z = [0001]);  same thermostat/barostat/rate as the composites.
variable        a0 equal 3.20
lattice         hcp ${a0}
region          box block 0 62 0 36 0 38 units lattice
create_box      3 box
include         ../potential.mod
create_atoms    1 box

thermo          1000
thermo_style    custom step temp pe press lx ly lz
fix             0 all box/relax aniso 0.0 vmax 0.001
min_style       cg
minimize        1.0e-12 1.0e-12 20000 200000
unfix           0

reset_timestep  0
timestep        0.001
velocity        all create 300.0 87287 dist gaussian
fix             1 all npt temp 300.0 300.0 0.1 aniso 0.0 0.0 5.0
run             20000
unfix           1
write_data      hcp_relaxed_300K.data

variable        temp equal 300.0
variable        srate equal 1.0e9
variable        erate equal v_srate/1.0e12
fix             2 all npt temp ${temp} ${temp} 0.1 __LAT1__ 0.0 0.0 5.0 __LAT2__ 0.0 0.0 5.0
fix             3 all deform 1 __AXIS__ erate ${erate} units box remap x
variable        L0 equal $(l__AXIS__)
variable        strain equal "(l__AXIS__ - v_L0)/v_L0"
variable        stress_x equal "-pxx/10000"
variable        stress_y equal "-pyy/10000"
variable        stress_z equal "-pzz/10000"
thermo_style    custom step v_strain v_stress_x v_stress_y v_stress_z temp pe press
fix             4 all print 100 "${strain} ${stress_x} ${stress_y} ${stress_z}" file tension_stress_strain.txt screen no title "Strain Stress_X_GPa Stress_Y_GPa Stress_Z_GPa"
# no per-atom dump: only the stress-strain curve is needed for the crystal reference
run             150000
"""
for ax in "xz":
    lat = ("y", "z") if ax == "x" else ("x", "y")
    t = HCP.replace("__AXIS__", ax).replace("__LAT1__", lat[0]).replace("__LAT2__", lat[1])
    add("hcp", f"hcp_{ax}", ["hcp"], dict(AXIS=ax), {"in.hcp_tension": t},
        dict(group=f"hcp_{ax}", xZn=0, xCa=0, vf=1.0, radius=0, axis=ax, purpose="perfect-crystal reference"))

# --------------------------------------------------------------------------------------------------
# Block 5  temp : temperature sweep, baseline composition (run_1 relaxed composite), 50 ps re-equilibration at each T
# --------------------------------------------------------------------------------------------------
for T in (250, 300, 350, 400):      # 300 K is a same-protocol control (50 ps re-equilibration) so the sweep is internally consistent
    add("temp", f"temp_{T}", ["copyrelaxed", "tension"],
        dict(SOURCE_RELAXED="../run_1/AN_composite_relaxed_300K.data", AXIS="x"),
        {"week4_tension.in": tension(temp=float(T), equil=50000)},
        dict(group=f"run1_T{T}", xZn=0.35, xCa=0.05, vf=0.007, radius=25.0, axis="x",
             purpose="temperature sweep (T=%d K)" % T))

# --------------------------------------------------------------------------------------------------
# Block 6  rate : 2e9 /s check at two further Zn levels (existing check was Composition 6, Ca-varying only)
# --------------------------------------------------------------------------------------------------
for run in (4, 3):
    add("rate", f"rate2e9_run{run}", ["copyrelaxed", "tension"],
        dict(SOURCE_RELAXED=f"../run_{run}/AN_composite_relaxed_300K.data", AXIS="x"),
        {"week4_tension.in": tension(srate=2.0e9)},
        dict(group=f"run{run}_2e9", xZn=COMP[run][0], xCa=COMP[run][1], vf=0.007, radius=COMP[run][2], axis="x",
             purpose="rate sensitivity at another Zn level"))

# --------------------------------------------------------------------------------------------------
# Block 7  clean : Compositions 1 and 9 repeated from a single uninterrupted quench (same seeds as the originals)
# --------------------------------------------------------------------------------------------------
for run in (1, 9):
    xzn, xca, rad = COMP[run]
    add("clean", f"clean_r{run}", ["quench", "build", "relax", "tension"], dict(AXIS="x"),
        {"week2_melt_quench.in": quench(xzn, xca),
         "build_an_composite.py": builder("quenched_matrix_300K.data", rad, None, 1000 + run),
         "in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
        dict(group=f"clean_r{run}", xZn=xzn, xCa=xca, vf=0.007, radius=rad, axis="x",
             purpose="uninterrupted-quench repeat of Composition %d" % run))

# --------------------------------------------------------------------------------------------------
# Block 8  indep_q : independent melt-quench realizations (new composition-assignment AND velocity seeds), then 12 A / 10 %.
#          These are the replicates that include glass-structure scatter.  Full pipeline = the expensive path.
# --------------------------------------------------------------------------------------------------
SEEDS = {2: (20123, 20321, 9902), 3: (30123, 30321, 9903)}
for rep, sd in SEEDS.items():
    for run in (5, 1, 6):
        xzn, xca, _ = COMP[run]
        add("indep_q", f"q{rep}_r{run}", ["quench", "build", "relax", "tension"], dict(AXIS="x"),
            {"week2_melt_quench.in": quench(xzn, xca, sd),
             "build_an_composite.py": builder("quenched_matrix_300K.data", VF_R, 0.10, 6000 + rep*10 + run),
             "in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
            dict(group=f"r{run}_vf10", xZn=xzn, xCa=xca, vf=0.10, radius=VF_R, axis="x",
                 purpose=f"independent quench realization {rep}"))

# --------------------------------------------------------------------------------------------------
# Block 9  large : 2x2x2 replica of the run_1 glass (~3.1 M atoms, ~41 nm), 12 A seeds at 5 and 15 %.  OPTIONAL, expensive.
# --------------------------------------------------------------------------------------------------
REP = """include         ../common.mod
read_data       ../run_1/quenched_matrix_300K.data
include         ../potential.mod
replicate       2 2 2
# fresh velocities + 50 ps NPT so thermal motion is not 8 identical copies
timestep        0.001
velocity        all create 300.0 31415 dist gaussian
fix             1 all npt temp 300.0 300.0 0.1 iso 0.0 0.0 5.0
thermo          1000
run             50000
write_data      quenched_matrix_300K.data
"""
for k, vf in enumerate([0.05, 0.15]):
    add("large", f"large_r1_vf{int(vf*100)}", ["replicate", "build", "relax", "tension"],
        dict(SOURCE_GLASS="../run_1/quenched_matrix_300K.data", AXIS="x"),
        {"in.replicate": REP, "build_an_composite.py": builder("quenched_matrix_300K.data", VF_R, vf, 7001 + k),
         "in.relax_composite": RELAX_TEMPLATE, "week4_tension.in": tension()},
        dict(group=f"large_vf{int(vf*100)}", xZn=0.35, xCa=0.05, vf=vf, radius=VF_R, axis="x",
             purpose="2x2x2 cell (shear-band size test)"), ntasks=160)

# --------------------------------------------------------------------------------------------------
os.makedirs("cases", exist_ok=True); os.makedirs("logs", exist_ok=True)
for b, d in blocks.items():
    open(f"cases/{b}.txt", "w").write("\n".join(d["cases"]) + "\n")
    open(f"cases/{b}.ntasks", "w").write(str(d["ntasks"]) + "\n")
with open("cases_manifest.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["case", "block", "group", "xZn", "xCa", "vf", "radius", "axis", "purpose"])
    w.writeheader(); w.writerows(manifest)
print(f"{len(manifest)} cases in {len(blocks)} blocks")
for b, d in blocks.items(): print(f"  {b:11s} {len(d['cases']):2d} cases  ntasks={d['ntasks']:3d}  {' '.join(d['cases'])}")
