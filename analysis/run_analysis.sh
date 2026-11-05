#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --output=analysis_%j.out
#SBATCH --error=analysis_%j.err
#SBATCH --partition=testq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#
# Generic wrapper for any standalone Python analysis script (NOT
# LAMMPS -- those already use their own sbatch scripts with mpirun).
# Runs the exact command given as arguments, on a compute node, never
# on the login/master node.
#
# Usage:
#   sbatch run_analysis.sh python3 analyze_d2min_trajectory.py vf_2/tension_trajectory.dump --erate_per_ps 0.001 --mode incremental
#
# testq's 1-hour limit comfortably covers every analysis script run
# so far (D2min: ~15-20 min per trajectory; extract_props.py,
# detect_shear_bands.py: well under a minute). If a future script
# genuinely needs longer, change --partition to mediumq and raise
# --time accordingly -- but default to testq for anything in this
# category.

set -euo pipefail
source ~/shap_env/bin/activate
echo "=== Running: $@ ==="
"$@"
echo "=== Complete ==="
