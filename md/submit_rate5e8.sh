#!/bin/bash
#SBATCH --job-name=an_rate5e8
#SBATCH --output=array_%j.out
#SBATCH --error=array_%j.err
#SBATCH --partition=longq
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=24-00:00:00
#
# GP round-2 active-learning point: strain_rate=5e8 /s (half baseline),
# otherwise identical to Composition 1. Half the rate means DOUBLE
# the steps (300,000 vs baseline's 150,000) are needed to reach the
# same final strain -- estimate ~24-28 hours, roughly 2x the T=310K/
# 2e9 cases' cost, not comparable to those in cheapness.
#
# Same reused-structure pattern: skips Stage 1-3 entirely, only
# Stage 4 (tension, modified rate/step count) runs.
#
# Requires: run_1/AN_composite_relaxed_300K.data and
# week4_tension_5e8.in, both in place before this can run.

set -euo pipefail

module purge
module load apps/anaconda3/2024.06
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/11Dec2020
export PYTHONNOUSERSITE=1

LMP_EXEC="/gpfs-home/p230674mt/lammps_src/src/lmp_mpi"

if [ ! -f "../run_1/AN_composite_relaxed_300K.data" ]; then
    echo "ERROR: ../run_1/AN_composite_relaxed_300K.data not found." >&2
    exit 1
fi
if [ ! -f "week4_tension_5e8.in" ]; then
    echo "ERROR: week4_tension_5e8.in not found in this directory." >&2
    exit 1
fi

cp ../run_1/AN_composite_relaxed_300K.data .

echo "=== Strain-rate=5e8/s Active-Learning Point (source: run_1) ==="
mpirun -np $SLURM_NTASKS $LMP_EXEC -in week4_tension_5e8.in -log log.4_tension_5e8

echo "=== Complete: rate5e8_1 ==="
