#!/bin/bash
#SBATCH -p longq
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J mix_${1}_${4}
cd /gpfs-home/p230674mt/W2_new/new_md_package/checks
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
LMP=/gpfs-home/p230674mt/lammps_src/src/lmp_mpi
# args: TAG XZN XCA SEED [MELT]
mpirun -np 8 $LMP -in in.mixing_check -var TAG $1 -var XZN $2 -var XCA $3 -var SEED $4 -var MELT ${5:-1} -screen none
