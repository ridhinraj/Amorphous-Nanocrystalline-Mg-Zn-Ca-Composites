#!/bin/bash
# run_case.sh <block> <line-number>   -- executed inside the SLURM job (or by hand: ./run_case.sh ca_matched 1)
# Runs the stages listed in <case>/case.cfg, skipping any stage whose output already exists, so a failed or
# timed-out job can simply be resubmitted.
set -euo pipefail
BLOCK="$1"; IDX="${2:-${SLURM_ARRAY_TASK_ID:-1}}"
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$ROOT"
source ./cluster.env
CASE=$(sed -n "${IDX}p" "cases/${BLOCK}.txt")
[ -n "$CASE" ] || { echo "no case on line $IDX of cases/${BLOCK}.txt" >&2; exit 1; }
NTASKS="${SLURM_NTASKS:-$(cat cases/${BLOCK}.ntasks)}"
if [ "${DRY_RUN:-0}" != "1" ]; then
  module purge; for m in $MODULES; do module load "$m"; done
  export PYTHONNOUSERSITE=1
  export PYTHONPATH="$PYMATGEN_LOCAL:${PYTHONPATH:-}"
fi
run() { if [ "${DRY_RUN:-0}" = "1" ]; then echo "[dry] $*"; else eval "$@"; fi; }
lmp() { run "$MPIRUN -np $NTASKS $LMP_EXEC -in $1 -log $2"; }

cd "$CASE"; source ./case.cfg
echo "=== $CASE  (block $BLOCK)  stages: $STAGES  ntasks=$NTASKS  $(date) ==="
tension_done() { [ -f log.4_tension ] && grep -q "Total wall time" log.4_tension; }

for S in $STAGES; do
  case $S in
    quench)      [ -f quenched_matrix_300K.data ] && { echo "skip quench"; continue; }
                 lmp week2_melt_quench.in log.1_melt_quench ;;
    replicate)   [ -f quenched_matrix_300K.data ] && { echo "skip replicate"; continue; }
                 [ -f "$SOURCE_GLASS" ] || { echo "missing $SOURCE_GLASS" >&2; exit 1; }
                 lmp in.replicate log.1_replicate ;;
    build)       [ -f MgZnCa_Composite_Ready.data ] && { echo "skip build"; continue; }
                 [ -z "${SOURCE_GLASS:-}" ] || [ -f "$SOURCE_GLASS" ] || { echo "missing $SOURCE_GLASS" >&2; exit 1; }
                 run "env -u LD_LIBRARY_PATH python3 build_an_composite.py 2>&1 | tee log.2_build" ;;
    copyglass)   [ -f MgZnCa_Composite_Ready.data ] && { echo "skip copyglass"; continue; }
                 [ -f "$SOURCE_GLASS" ] || { echo "missing $SOURCE_GLASS" >&2; exit 1; }
                 run "cp $SOURCE_GLASS MgZnCa_Composite_Ready.data" ;;
    relax)       [ -f AN_composite_relaxed_300K.data ] && { echo "skip relax"; continue; }
                 lmp in.relax_composite log.3_relax ;;
    copyrelaxed) [ -f AN_composite_relaxed_300K.data ] && { echo "skip copyrelaxed"; continue; }
                 [ -f "$SOURCE_RELAXED" ] || { echo "missing $SOURCE_RELAXED" >&2; exit 1; }
                 run "cp $SOURCE_RELAXED AN_composite_relaxed_300K.data" ;;
    tension)     tension_done && { echo "skip tension (finished)"; continue; }
                 lmp week4_tension.in log.4_tension ;;
    hcp)         tension_done && { echo "skip hcp (finished)"; continue; }
                 lmp in.hcp_tension log.4_tension ;;
    *) echo "unknown stage $S" >&2; exit 1 ;;
  esac
done
echo "=== done $CASE $(date) ==="
