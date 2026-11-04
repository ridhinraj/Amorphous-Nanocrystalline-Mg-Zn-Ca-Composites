#!/bin/bash
# verify_inputs.sh -- run on CLUSTER 1 in the project root after unpacking the bundle and running setup_new_md.py.
source ./cluster.env 2>/dev/null || { echo "cluster.env missing"; exit 1; }
bad=0; ok(){ printf "  ok    %s\n" "$1"; }; no(){ printf "  FAIL  %s\n" "$1"; bad=1; }; warn(){ printf "  WARN  %s\n" "$1"; }
echo "[1] shared inputs"
for f in common.mod potential.mod potentials/library.meam potentials/CaZnMg.meam md_templates.py setup_new_md.py run_case.sh submit_block.sh; do
  [ -f "$f" ] && ok "$f" || no "$f missing"; done
echo "[2] source structures (atom count / box edge from the data-file header)"
for f in run_1/quenched_matrix_300K.data run_5/quenched_matrix_300K.data run_6/quenched_matrix_300K.data \
         run_1/AN_composite_relaxed_300K.data run_3/AN_composite_relaxed_300K.data run_4/AN_composite_relaxed_300K.data; do
  if [ -f "$f" ]; then
    n=$(awk 'NR<20 && /atoms/{print $1; exit}' "$f"); lx=$(awk 'NR<20 && /xlo xhi/{printf "%.2f",$2-$1; exit}' "$f")
    case "$f" in *quenched*) [ "$n" = "389344" ] && ok "$f  $n atoms, lx=${lx} A" || warn "$f  $n atoms (expected 389344), lx=${lx} A";;
                 *)          ok "$f  $n atoms, lx=${lx} A";; esac
  else no "$f missing"; fi; done
echo "[3] LAMMPS on this cluster"
if [ -x "$LMP_EXEC" ]; then
  v=$("$LMP_EXEC" -h 2>&1 | grep -m1 -E "^Large-scale|LAMMPS \(" ); echo "  $v"
  "$LMP_EXEC" -h 2>&1 | grep -q -w meam && ok "pair_style meam present" || no "pair_style meam NOT built (needs the MEAM package)"
  "$LMP_EXEC" -h 2>&1 | grep -q -w zbl  && ok "pair_style zbl present"  || no "pair_style zbl NOT built"
else no "LMP_EXEC not executable: $LMP_EXEC  (edit cluster.env)"; fi
echo "  note: if a run stops with 'Fix halt command before simulation box is defined', run ./patch_for_new_lammps.sh once"
echo "[4] python for the build stage"
for m in $MODULES; do module load "$m" 2>/dev/null; done
PYTHONNOUSERSITE=1 PYTHONPATH="$PYMATGEN_LOCAL:${PYTHONPATH:-}" env -u LD_LIBRARY_PATH python3 -c "import numpy, scipy, pandas; from pymatgen.io.lammps.data import LammpsData; print('  ok    numpy/scipy/pandas/pymatgen import')" 2>&1 | tail -1
echo "[5] scheduler / generated cases"
command -v sbatch >/dev/null && ok sbatch || no "sbatch not found"
[ -f cases_manifest.csv ] && ok "cases generated: $(($(wc -l < cases_manifest.csv)-1)) cases" || no "run python3 setup_new_md.py first"
[ $bad -eq 0 ] && echo "ALL CHECKS PASSED" || echo "FIX THE FAIL LINES ABOVE"
exit $bad
