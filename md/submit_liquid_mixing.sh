#!/bin/bash
#SBATCH --job-name=phase1_liquid_mix
#SBATCH --output=phase1_validation/liquid_%a_%A.out
#SBATCH --error=phase1_validation/liquid_%a_%A.err
#SBATCH --partition=testq
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --array=1-4
#
# The 4 compositions (Mg, Zn, Ca, MgZnCa) are fully independent of
# each other and of lattice_relax.in/gsfe_basal.in -- runs as its
# own array, in parallel with submit_lattice_gsfe.sh if you like.
#
# All 4 scripts stay flat inside phase1_validation/ (no
# subdirectories per array task) because potential.mod's relative
# path to ../potentials/ depends on the working directory being
# exactly one level below project_root -- same constraint as
# run_1..run_10 and vf_1..vf_3.

set -euo pipefail

module purge
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/11Dec2020

LMP_EXEC="/gpfs-home/p230674mt/lammps_src/src/lmp_mpi"

cd "$SLURM_SUBMIT_DIR/phase1_validation"

case $SLURM_ARRAY_TASK_ID in
    1) SCRIPT=liquid_Mg.in ;;
    2) SCRIPT=liquid_Zn.in ;;
    3) SCRIPT=liquid_Ca.in ;;
    4) SCRIPT=liquid_MgZnCa.in ;;
    *) echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== Liquid mixing: $SCRIPT (array task $SLURM_ARRAY_TASK_ID) ==="
mpirun -np $SLURM_NTASKS $LMP_EXEC -in "$SCRIPT" -log "log.$SCRIPT"
grep PHASE1_RESULT "log.$SCRIPT"
echo "=== Done: $SCRIPT ==="
