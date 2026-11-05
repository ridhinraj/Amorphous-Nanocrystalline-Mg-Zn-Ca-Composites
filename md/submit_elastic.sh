#!/bin/bash
#SBATCH --job-name=an_elastic
#SBATCH --output=array_%j.out
#SBATCH --error=array_%j.err
#SBATCH --partition=testq
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#
# Elastic constants calculation, adapted from LAMMPS's standard
# ELASTIC example. Given the small system size and short runtime of
# every other Phase 1 validation calculation in this project
# (~20 seconds each on testq), this should be similarly fast --
# testq's 1-hour limit gives ample margin.
#
# Run from within phase1_validation/ (or wherever this project's
# ../common.mod and ../potential.mod are one level up).

set -euo pipefail

module purge
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/11Dec2020

LMP_EXEC="/gpfs-home/p230674mt/lammps_src/src/lmp_mpi"

mpirun -np $SLURM_NTASKS $LMP_EXEC -in in.elastic -log log.elastic

echo "=== Elastic constants calculation complete ==="
echo "Check log.elastic for the SELF-CHECK lines before trusting results:"
echo "  1. All diagonal Cii (C11,C22,C33,C44,C55,C66) should be positive"
echo "  2. C66 should approximately equal (C11-C12)/2 (hexagonal symmetry)"
echo "  3. Matrix should be approximately symmetric (Cij ~ Cji)"
echo "  4. Compare against NIST reference: C11~63, C33~70, C12~25, C13~21, C44~17, C66~19 GPa"
