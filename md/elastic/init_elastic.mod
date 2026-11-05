# ============================================================
# init_elastic.mod
# Setup for elastic constants calculation, adapted from LAMMPS's
# standard ELASTIC example (examples/ELASTIC in the LAMMPS
# distribution) for this project's 2NN MEAM + ZBL potential.
#
# IMPORTANT: uses "lattice custom" with our own CONFIRMED relaxed
# c/a = 1.619 (from Phase 1 validation), not LAMMPS's built-in
# "lattice hcp" keyword, which hardcodes the ideal c/a = sqrt(8/3)
# = 1.633. Using the built-in keyword would silently build a
# structure inconsistent with this potential's actual relaxed
# geometry. The basis vectors below were numerically verified
# (independently, in Python) to reproduce Mg's known experimental
# density (1.742 g/cm^3 computed vs ~1.74 g/cm^3 real) and the
# correct near-a nearest-neighbor distance before being used here.
# ============================================================

variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5

units metal
variable cfac equal 1.0e-4
variable cunits string GPa

variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2

boundary p p p
atom_style atomic
atom_modify map array

# HCP Mg, orthogonal 4-atom cell, confirmed a=3.209 Ang, c/a=1.619
variable a0 equal 3.209
variable ca_ratio equal 1.619

lattice custom ${a0} &
    a1 1.0 0.0 0.0 &
    a2 0.0 1.7320508 0.0 &
    a3 0.0 0.0 ${ca_ratio} &
    basis 0.0 0.0 0.0 &
    basis 0.5 0.5 0.0 &
    basis 0.5 0.8333333333 0.5 &
    basis 0.0 0.3333333333 0.5

region box block 0 6 0 6 0 6 units lattice
create_box 3 box
create_atoms 1 box

mass 1 24.305
mass 2 65.38
mass 3 40.078

include ../potential.mod

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes

# Small random jiggle to avoid atoms sitting exactly on a saddle
# point of the potential energy surface (standard practice in the
# original ELASTIC example, not specific to this adaptation).
displace_atoms all random ${atomjiggle} ${atomjiggle} ${atomjiggle} 87287 units box

# Relax the full cell (shape + volume) before applying any strain,
# so the reference state is the true equilibrium structure -- this
# is where we expect to recover a=3.209, c/a=1.619 as a
# self-consistency check (see in.elastic for the explicit check).
fix 3 all box/relax aniso 0.0 vmax 0.001
min_style cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix 3

# Store reference box dimensions for use in displace_elastic.mod
variable lx0 equal "lx"
variable ly0 equal "ly"
variable lz0 equal "lz"

# Store reference stress tensor (should be near zero after relaxation)
variable pxx0 equal "pxx"
variable pyy0 equal "pyy"
variable pzz0 equal "pzz"
variable pxy0 equal "pxy"
variable pxz0 equal "pxz"
variable pyz0 equal "pyz"

# Self-consistency check: print the actual relaxed lattice
# parameters so they can be directly compared against the confirmed
# Phase 1 values (a=3.209, c/a=1.619) before trusting anything
# downstream. If these don't match closely, STOP and investigate
# before proceeding to the strain calculations.
print "SELF-CHECK: relaxed lx = ${lx0} (expect close to 6*3.209 = 19.254)"
print "SELF-CHECK: relaxed lz/6 / a0 -> c/a should be close to 1.619"
variable check_a equal "lx/6"
variable check_ca equal "(lz/6)/v_check_a"
print "SELF-CHECK: computed a = ${check_a} Ang"
print "SELF-CHECK: computed c/a = ${check_ca}"

write_restart restart.equil
