#!/bin/bash
#SBATCH --job-name=A4_glass0
#SBATCH --output=A4_%x_%j.out
#SBATCH --error=A4_%x_%j.err
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#
# A4: pure-glass (no crystal) controls.  Answers two questions without any crystal present:
#   (a) incremental D2min -- how many simultaneously active rearranging clusters does the plain glass show (the "0 % VF" point of Table 9)?
#   (b) per-atom stress / Voronoi-volume analysis -- does the Ca-vs-Zn stress difference, the structure/volume explained fraction and
#       the Ca load-partition factor exist in a glass with no nanocrystal?
#
# Usage (from the folder that contains glass0_r5/ ... and the scripts; one job per case):
#   sbatch run_A4_glass0.sh glass0_r5 0.0703
#   sbatch run_A4_glass0.sh glass0_r1 0.0657
#   sbatch run_A4_glass0.sh glass0_r6 0.0722
# 2nd argument = peak strain from new_md_properties_final.csv.  Optional env: SCRIPTS=<folder holding the .py files and voronoi_volume/> (default: current folder).
# Runs on 1 core (+ threads); ~20 min D2min + ~10 min Voronoi (3 frames) + ~10 min regression.  Use --partition=shortq or mediumq as slots allow.
set -euo pipefail
CASE=${1:?case folder, e.g. glass0_r5}; PEAK=${2:?peak strain, e.g. 0.0703}
SCRIPTS=${SCRIPTS:-$PWD}
source ~/shap_env/bin/activate
for f in analyze_d2min_trajectory.py compute_d2min.py combined_regression.py fastlib.py multiframe_solute_shap.py voronoi_volume/select_frames.py voronoi_volume/voronoi_py.py voronoi_volume/voronoi_merge_v2.py; do
  [ -f "$SCRIPTS/$f" ] || { echo "MISSING $SCRIPTS/$f"; exit 1; }; done
[ -f "$CASE/tension_trajectory.dump" ] || { echo "MISSING $CASE/tension_trajectory.dump"; exit 1; }
export PYTHONPATH="$SCRIPTS:${PYTHONPATH:-}"
cd "$CASE"
echo "=== [1/4] incremental D2min  ($(date)) ==="
python3 "$SCRIPTS/analyze_d2min_trajectory.py" tension_trajectory.dump --erate_per_ps 0.001 --mode incremental | tee A4_d2min_incremental.log
echo "=== [2/4] frames 0, peak ($PEAK), 0.15 ==="
python3 "$SCRIPTS/voronoi_volume/select_frames.py" tension_trajectory.dump --list | head -5
python3 "$SCRIPTS/voronoi_volume/select_frames.py" tension_trajectory.dump --strains 0 "$PEAK" 0.15 \
        --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 -o frames_selected.dump
echo "=== [3/4] Voronoi volumes  ($(date)) ==="
python3 "$SCRIPTS/voronoi_volume/voronoi_py.py" frames_selected.dump -o volumes.dump | tee A4_voronoi.log
python3 "$SCRIPTS/voronoi_volume/voronoi_merge_v2.py" frames_selected.dump volumes.dump --label "$CASE" \
        --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 --stress-file tension_stress_strain.txt --save-atoms "${CASE}_atoms.npz" | tee A4_merge.log
echo "=== [4/4] regression with volume-normalised stress  ($(date)) ==="
for S in 0 "$PEAK" 0.15; do
  python3 "$SCRIPTS/combined_regression.py" tension_trajectory.dump "$S" --erate_per_ps 0.001 --volumes volumes.dump --out "A4_cr_${CASE}_${S}.csv" | tee "A4_cr_${CASE}_${S}.log"
done
echo "=== DONE $CASE ($(date)) ==="
echo "Send back: ${CASE}/A4_*.log  ${CASE}/A4_cr_*.csv  ${CASE}/*_d2min_incremental_summary.csv  ${CASE}/${CASE}_species_summary.csv (if written)"
