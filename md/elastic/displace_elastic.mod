# ============================================================
# displace_elastic.mod
# Applies ONE small strain (magnitude ${up}, sign ${updown}: +1 or
# -1) in direction ${dir} (1-6 = xx,yy,zz,yz,xz,xy in Voigt
# notation), re-minimizes atomic positions at the strained cell
# shape, and records the resulting stress tensor. Called twice per
# direction (once with updown=+1, once with updown=-1) by
# in.elastic, so that Cij can be computed via central difference
# rather than a one-sided (less accurate) difference.
#
# Included from a FRESH read of the equilibrium restart each time,
# so each of the 6 directions starts from the identical, already-
# validated equilibrium structure (not from a previous direction's
# strained state).
# ============================================================

clear
box tilt large
read_restart restart.equil
include ../potential.mod

# The box was built as a plain orthogonal box (via a "block" region
# in init_elastic.mod), so it cannot have shear tilt factors applied
# directly -- LAMMPS confirmed this explicitly on first real test:
# "ERROR: Cannot change box tilt factors for orthogonal box" on the
# first shear direction (yz, dir=4) after the three normal-strain
# directions (xx/yy/zz, dir=1-3) succeeded. Converting to triclinic
# here, unconditionally, is harmless for the normal-strain directions
# too (adds zero-valued tilt factors initially) and is required
# before any yz/xz/xy change_box call can succeed.
change_box all triclinic

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes

# Defensive: explicitly delete before redefining, since this script
# is included 12 times total (6 directions x 2 signs) and each of
# these variable names gets reused every time. LAMMPS documentation
# states equal-style variables CAN be redefined directly, but since
# this could not be verified by actually running LAMMPS in the
# environment this was written in, deleting first removes any doubt
# at zero cost (delete is a no-op if the variable doesn't exist yet).
variable delta delete
variable deltaxy delete
variable deltaxz delete
variable deltayz delete
variable delta equal ${up}*${updown}*${len0}
variable deltaxy equal ${up}*${updown}*ly
variable deltaxz equal ${up}*${updown}*lz
variable deltayz equal ${up}*${updown}*lz

if "${dir} == 1" then &
    "change_box all x delta 0 ${delta} remap units box"
if "${dir} == 2" then &
    "change_box all y delta 0 ${delta} remap units box"
if "${dir} == 3" then &
    "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
    "change_box all yz delta ${deltayz} remap units box"
if "${dir} == 5" then &
    "change_box all xz delta ${deltaxz} remap units box"
if "${dir} == 6" then &
    "change_box all xy delta ${deltaxy} remap units box"

minimize ${etol} ${ftol} ${maxiter} ${maxeval}

variable tmp delete
variable pxx1 delete
variable pyy1 delete
variable pzz1 delete
variable pyz1 delete
variable pxz1 delete
variable pxy1 delete
variable tmp equal pxx
variable pxx1 equal ${tmp}
variable tmp delete
variable tmp equal pyy
variable pyy1 equal ${tmp}
variable tmp delete
variable tmp equal pzz
variable pzz1 equal ${tmp}
variable tmp delete
variable tmp equal pyz
variable pyz1 equal ${tmp}
variable tmp delete
variable tmp equal pxz
variable pxz1 equal ${tmp}
variable tmp delete
variable tmp equal pxy
variable pxy1 equal ${tmp}
