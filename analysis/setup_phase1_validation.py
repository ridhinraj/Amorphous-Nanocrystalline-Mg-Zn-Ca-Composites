import os

# ============================================================
# setup_phase1_validation.py (corrected)
#
# Fixes two issues found from actually running this on testq:
#
# 1. BUG: the liquid-mixing "production" phase had no active
#    integrator. fix 1 (npt) was unfixed before the production
#    run, and fix ave/time was left to average a value over a
#    run with nothing moving the atoms -- LAMMPS printed
#    "WARNING: No fixes with time integration, atoms won't move"
#    and every sampled point in the output file was identical,
#    i.e. one frozen snapshot repeated, not a real average.
#    Fixed by keeping a second npt fix active through production.
#
# 2. REALITY CHECK: testq granted 4 MPI tasks for a 16-task
#    request (1x2x2 grid in the log), not 16 -- a 4x slowdown
#    from what the original sizing assumed. Combined with the
#    bug above, task 1 (Mg) took 42m41s against a 1-hour limit,
#    with zero restart capability if it had been killed.
#    Step counts below are cut to 15000+5000 (from 30000+10000)
#    to rebuild real safety margin under the confirmed 4-core
#    reality, not the assumed 16-core one.
# ============================================================

os.makedirs("phase1_validation", exist_ok=True)

# ------------------------------------------------------------
# 1. Lattice relaxation -- pure Mg, HCP (unchanged, already
#    confirmed working: a=3.209 Ang, c/a=1.619, Ecoh=-1.549 eV/atom)
# ------------------------------------------------------------
lattice_relax = """include         ../common.mod

variable        a0 equal 3.19          # starting guess only -- box/relax finds the true value
lattice         hcp ${a0}
region          box block 0 10 0 10 0 10 units lattice
create_box      3 box
create_atoms    1 box                  # type 1 = Mg only

include         ../potential.mod

compute         peatom all pe/atom

fix             1 all box/relax aniso 0.0 vmax 0.001
min_style       cg
minimize        1.0e-14 1.0e-14 20000 200000
unfix           1

variable        a_relaxed equal "lx/10"
variable        c_relaxed equal "lz/10"
variable        Ecoh      equal "pe/atoms"

print           "==================================================="
print           "Relaxed Mg lattice:  a = ${a_relaxed} Ang   c = ${c_relaxed} Ang"
print           "c/a ratio          = $(v_c_relaxed/v_a_relaxed)"
print           "Cohesive energy     = ${Ecoh} eV/atom"
print           "==================================================="
print           "Compare against literature: a ~ 3.19 Ang, c/a ~ 1.624, Ecoh ~ -1.51 eV/atom"
print           "PHASE1_RESULT a_relaxed=${a_relaxed} c_relaxed=${c_relaxed} Ecoh=${Ecoh}"

write_data      Mg_hcp_relaxed.data
write_restart   Mg_hcp_relaxed.restart
"""
with open("phase1_validation/lattice_relax.in", "w") as f:
    f.write(lattice_relax)

# ------------------------------------------------------------
# 2. Basal GSFE -- unchanged from the confirmed-working version.
#    UPDATE bshift to 3.20906499038793 (the real a_relaxed from
#    the run already completed) before this next run, replacing
#    the 3.19 placeholder used previously.
# ------------------------------------------------------------
gsfe_basal = """# ============================================================
# Basal-plane GSFE, pure Mg. Requires lattice_relax.in to have
# already run in this same directory (reads its restart file).
#
# bshift below is the ACTUAL relaxed a0 from lattice_relax.in's
# PHASE1_RESULT (3.20906499038793), not the earlier 3.19 guess.
# ============================================================

variable        bshift  equal 3.20906499038793
variable        nsteps  equal 20

log             gsfe.log

label           gsfe_loop
variable        n loop 0 ${nsteps}

  clear
  include       ../common.mod
  read_restart  Mg_hcp_relaxed.restart
  include       ../potential.mod

  change_box    all boundary p p f
  change_box    all z final $(zlo-15) $(zhi+15) units box

  variable      zmid equal "(zlo+zhi)/2.0"
  region        lower block EDGE EDGE EDGE EDGE EDGE ${zmid} units box
  region        upper block EDGE EDGE EDGE EDGE ${zmid} EDGE units box
  group         lower region lower
  group         upper region upper

  variable      dx equal "v_bshift*v_n/v_nsteps"
  displace_atoms upper move ${dx} 0.0 0.0 units box

  fix           1 lower setforce 0.0 0.0 NULL
  fix           2 upper setforce 0.0 0.0 NULL
  min_style     cg
  minimize      1.0e-14 1.0e-14 10000 100000

  variable      etot equal pe
  variable      area equal "lx*ly"
  print         "GSFE_POINT n=${n} frac=$(v_n/v_nsteps) E=${etot} Area=${area}"

next            n
jump            SELF gsfe_loop

# Post-process: gamma(n) = (E(n) - E(0)) / Area, in eV/Ang^2;
# multiply by 16021.77 to convert to mJ/m^2. The confirmed result
# from the last run: peak (n=10, frac=0.5) gamma_us ~ 321 mJ/m^2,
# vs. ~95-120 mJ/m^2 from DFT literature for Mg basal <a> slip --
# a ~3x gap worth checking against the Jang/Kim/Lee (2018) Mg-Zn
# paper's own published basal GSFE figure before trusting further.
"""
with open("phase1_validation/gsfe_basal.in", "w") as f:
    f.write(gsfe_basal)

# ------------------------------------------------------------
# 3. Liquid mixing enthalpy -- CORRECTED: fix 2 (npt) now stays
#    active through the production/averaging run instead of
#    leaving the system with no integrator. Step counts cut to
#    15000+5000 (from 30000+10000) for real safety margin at the
#    confirmed 4 cores, not the assumed 16.
# ------------------------------------------------------------
LIQUID_TEMPLATE = """include         ../common.mod

# NOTE: Ttarget below is a placeholder (1000 K), not confirmed
# against the specific temperature used in the original Jang,
# Seol & Lee (2019) CALPHAD liquid-mixing assessment -- that
# number is buried in a figure in the paper and wasn't found via
# search. Check the paper directly and adjust before treating
# results as a real validation rather than a sanity check.
variable        Ttarget equal 1000.0
variable        xZn     equal __XZN__
variable        xCa     equal __XCA__
variable        tag     string "__TAG__"

lattice         fcc 4.5
region          box block 0 12 0 12 0 12 units lattice   # 6912 atoms
create_box      3 box
create_atoms    1 box

include         ../potential.mod

if "${xZn} > 0.0" then &
  "set type 1 type/fraction 2 ${xZn} 12345"
if "${xCa} > 0.0" then &
  "set type 1 type/fraction 3 $(v_xCa/(1.0-v_xZn)) 54321"

velocity        all create ${Ttarget} 9834 dist gaussian

fix             1 all npt temp ${Ttarget} ${Ttarget} 0.1 iso 0.0 0.0 1.0
timestep        0.001
thermo          2000
thermo_style    custom step temp pe press vol density
run             15000   # melt + equilibrate -- confirmed-safe length at 4 cores; grow only if you confirm more cores are actually granted
unfix           1

# Production/averaging run -- fix 2 keeps NPT integration active
# through this whole window. Without an active integrator here,
# fix ave/time just re-samples one frozen configuration repeatedly
# (confirmed from the previous run's log: identical thermo values
# at every sampled step, plus the "No fixes with time integration"
# warning).
reset_timestep  0
fix             2 all npt temp ${Ttarget} ${Ttarget} 0.1 iso 0.0 0.0 1.0
variable        peatom equal "pe/atoms"
fix             avg all ave/time 100 25 2500 v_peatom file liquid_pe_${tag}.dat
run             5000    # production average -- confirmed-safe length at 4 cores
unfix           2

print           "PHASE1_RESULT tag=${tag} xMg=$(1.0-v_xZn-v_xCa) xZn=${xZn} xCa=${xCa}"
print           "Mean PE/atom written to liquid_pe_${tag}.dat"
"""

liquid_compositions = [
    # tag,      xZn,   xCa
    ("Mg",      0.0,   0.0),
    ("Zn",      1.0,   0.0),
    ("Ca",      0.0,   1.0),
    ("MgZnCa",  0.35,  0.05),
]

for tag, xZn, xCa in liquid_compositions:
    content = (
        LIQUID_TEMPLATE
        .replace("__XZN__", str(xZn))
        .replace("__XCA__", str(xCa))
        .replace("__TAG__", tag)
    )
    with open(f"phase1_validation/liquid_{tag}.in", "w") as f:
        f.write(content)

# ------------------------------------------------------------
# 4. Elastic constants -- unchanged, still manual (see reasoning
#    in the file itself).
# ------------------------------------------------------------
elastic_readme = """# Elastic constants (Cij) -- manual step, not part of submit_testq.sh

Reuses LAMMPS's bundled ELASTIC example instead of custom code.

1. Locate it: <lammps_source>/examples/ELASTIC/ (or download the 4
   files -- init.mod, potential.mod, displace.mod, elastic.in -- from
   docs.lammps.org, search "ELASTIC").

2. Replace that example's potential.mod with this content:

   pair_style      hybrid/overlay meam zbl 1.0 2.0
   pair_coeff      * * meam ../potentials/library.meam Ca Zn Mg ../potentials/CaZnMg.meam Mg Zn Ca
   pair_coeff      1 1 zbl 12.0 12.0
   pair_coeff      1 2 zbl 12.0 30.0
   pair_coeff      1 3 zbl 12.0 20.0
   pair_coeff      2 2 zbl 30.0 30.0
   pair_coeff      2 3 zbl 30.0 20.0
   pair_coeff      3 3 zbl 20.0 20.0

   mass            1 24.305
   mass            2 65.38
   mass            3 40.078

3. In that example's init.mod, set:
     lattice hcp 3.20906499038793   (confirmed a_relaxed from lattice_relax.in)
     create_atoms 1 box

4. Run: lammps -in elastic.in
   Prints the 6x6 elastic stiffness matrix -- expect 5 independent
   nonzero constants for HCP Mg (C11, C12, C13, C33, C44).

Small enough to run in the same testq window as everything else here,
but kept manual and separate since it's a different (external,
pre-tested) script rather than one generated by this project.
"""
with open("phase1_validation/ELASTIC_CONSTANTS.md", "w") as f:
    f.write(elastic_readme)

print("Generated phase1_validation/: lattice_relax.in, gsfe_basal.in, "
      "liquid_Mg.in, liquid_Zn.in, liquid_Ca.in, liquid_MgZnCa.in, "
      "ELASTIC_CONSTANTS.md")
print("liquid_*.in: integrator bug fixed, steps cut to 15000+5000 for 4-core safety margin.")
print("gsfe_basal.in: bshift updated to the confirmed a_relaxed value.")
