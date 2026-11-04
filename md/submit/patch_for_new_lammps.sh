#!/bin/bash
# Only needed if cluster 1's LAMMPS is newer than the 11Dec2020 build used on cluster 2 and stops with
#   "ERROR: Fix halt command before simulation box is defined"
# Newer versions require the box to exist before 'fix halt'.  common.mod is included BEFORE the box exists in every
# script, potential.mod AFTER it, so the halt safety-net is moved from the first file to the second.  Nothing else changes.
set -euo pipefail
[ -f common.mod ] && [ -f potential.mod ] || { echo "run from the project root"; exit 1; }
grep -q "halt_check" common.mod || { echo "common.mod has no halt_check line (already patched?)"; exit 0; }
[ -f common.mod.orig ] || cp common.mod common.mod.orig; [ -f potential.mod.orig ] || cp potential.mod potential.mod.orig
grep -v -E "^variable current_temp|^fix halt_check" common.mod > common.mod.tmp && mv common.mod.tmp common.mod
{ echo; echo "# moved here from common.mod by patch_for_new_lammps.sh (box must exist first)";
  echo "variable current_temp equal temp"; echo "fix halt_check all halt 200 v_current_temp > 5000 error hard"; } >> potential.mod
echo "patched: originals kept as common.mod.orig / potential.mod.orig"
