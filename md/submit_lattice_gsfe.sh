#!/bin/bash
#SBATCH --job-name=phase1_lattice_gsfe
#SBATCH --output=phase1_validation/lattice_gsfe_%j.out
#SBATCH --error=phase1_validation/lattice_gsfe_%j.err
#SBATCH --partition=testq
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#
# Both scripts here are small (few thousand atoms, minimizations
# only, no long MD) and should finish quickly -- ntasks=16 is a
# conservative default given testq's actual node config is unknown
# to me; adjust to whatever testq typically offers.
#
# gsfe_basal.in reads Mg_hcp_relaxed.restart, so lattice_relax.in
# MUST run first -- that's why these are chained in one job rather
# than two, unlike the independent liquid-mixing compositions.

set -euo pipefail

module purge
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/11Dec2020

LMP_EXEC="/gpfs-home/p230674mt/lammps_src/src/lmp_mpi"

cd "$SLURM_SUBMIT_DIR/phase1_validation"

echo "=== 1. Lattice relaxation ==="
mpirun -np $SLURM_NTASKS $LMP_EXEC -in lattice_relax.in -log log.lattice_relax
grep PHASE1_RESULT log.lattice_relax

echo "=== 2. Basal GSFE ==="
mpirun -np $SLURM_NTASKS $LMP_EXEC -in gsfe_basal.in -log log.gsfe_basal
grep GSFE_POINT gsfe.log | head -3
grep GSFE_POINT gsfe.log | tail -3

echo "=== Phase 1 lattice/GSFE complete ==="
echo "Reminder: gsfe_basal.in used a placeholder bshift (3.19) -- check"
echo "the a_relaxed value above and update bshift in gsfe_basal.in if"
echo "it differs meaningfully, then rerun for a corrected GSFE curve."
