#!/bin/bash
# submit_chain.sh [--partition P] [--time D-HH:MM:SS] [--name N] [--dry] block:line [block:line ...]
#   ONE job that runs the listed cases one after another (all with the same MPI rank count).
set -euo pipefail
source ./cluster.env
NAME=chain; DRY=0; NODES=""; LIST=()
while [ $# -gt 0 ]; do
  case "$1" in
    --partition) PARTITION="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --name) NAME="$2"; shift 2;;
    --nodes) NODES="$2"; shift 2;;
    --dry) DRY=1; shift;;
    *:*) LIST+=("$1"); shift;;
    *) echo "unknown argument $1"; exit 1;;
  esac
done
[ ${#LIST[@]} -gt 0 ] || { echo "give at least one block:line"; exit 1; }
NT=""
for c in "${LIST[@]}"; do
  b=${c%%:*}; i=${c##*:}
  [ -f "cases/${b}.txt" ] || { echo "no cases/${b}.txt"; exit 1; }
  [ -n "$(sed -n "${i}p" "cases/${b}.txt")" ] || { echo "no line $i in cases/${b}.txt"; exit 1; }
  n=$(cat "cases/${b}.ntasks"); [ -z "$NT" ] && NT=$n
  [ "$n" = "$NT" ] || { echo "mixed rank counts ($b uses $n, first case uses $NT)"; exit 1; }
  echo "  $c -> $(sed -n "${i}p" "cases/${b}.txt")"
done
mkdir -p logs
WRAP="for c in ${LIST[*]}; do b=\${c%%:*}; i=\${c##*:}; ./run_case.sh \$b \$i || echo FAILED \$c; done"
CMD=(sbatch --job-name="$NAME" --partition="$PARTITION" --ntasks="$NT" ${NODES:+--nodes=$NODES} --cpus-per-task=1 --time="$TIME"
     --output="logs/${NAME}_%j.out" --error="logs/${NAME}_%j.err")
[ -n "${EXCLUDE:-}" ] && CMD+=(--exclude="$EXCLUDE")
CMD+=(--wrap="$WRAP")
if [ "$DRY" = 1 ]; then printf '%q ' "${CMD[@]}"; echo; else "${CMD[@]}"; fi
