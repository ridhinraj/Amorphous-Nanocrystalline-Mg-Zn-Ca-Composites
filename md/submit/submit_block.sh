#!/bin/bash
# submit_block.sh <block> [--range 1-3] [--partition mediumq] [--time 2-00:00:00] [--dry]
#   e.g.  ./submit_block.sh glass0                       all cases of the block, PARTITION from cluster.env
#         ./submit_block.sh pack_reps --range 1-3 --partition mediumq
#         ./submit_block.sh pack_reps --range 4-4 --partition longq     (one case on the long queue)
# One SLURM array task per case listed in cases/<block>.txt.  Settings come from cluster.env.
set -euo pipefail
BLOCK="${1:?usage: submit_block.sh <block> [--range a-b] [--partition P] [--time D-HH:MM:SS] [--dry]}"; shift || true
source ./cluster.env
N=$(wc -l < "cases/${BLOCK}.txt"); NT=$(cat "cases/${BLOCK}.ntasks")
RANGE="1-${N}"; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --range) RANGE="$2"; shift 2;;
    --partition) PARTITION="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --dry) DRY=1; shift;;
    *) echo "unknown option $1"; exit 1;;
  esac
done
CMD=(sbatch --job-name="md_${BLOCK}" --partition="$PARTITION" --ntasks="$NT" --cpus-per-task=1 --time="$TIME"
     --array="${RANGE}" --output="logs/${BLOCK}_%a_%A.out" --error="logs/${BLOCK}_%a_%A.err")
[ -n "${EXCLUDE:-}" ] && CMD+=(--exclude="$EXCLUDE")
CMD+=(--wrap="./run_case.sh ${BLOCK}")
mkdir -p logs
if [ "$DRY" = 1 ]; then echo "${CMD[@]}"; else "${CMD[@]}"; fi
