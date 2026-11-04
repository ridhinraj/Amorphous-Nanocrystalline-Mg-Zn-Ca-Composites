#!/bin/bash
#SBATCH -p mediumq
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J mixall
cd /gpfs-home/p230674mt/W2_new/new_md_package/checks
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
LMP=/gpfs-home/p230674mt/lammps_src/src/lmp_mpi
run() { mpirun -np 4 $LMP -in in.mixing_check -var TAG $1 -var XZN $2 -var XCA $3 -var SEED $4 -var MELT ${5:-1} -screen none; }
# seed 1 first (the answer), with the old-protocol Ca control
run Mg 0.0 0.0 1 & run Zn 1.0 0.0 1 & run Ca 0.0 1.0 1 & run MgZnCa 0.35 0.05 1 & wait
run Ca_oldprotocol 0.0 1.0 1 0 &
for SEED in 2 3; do
  run Mg 0.0 0.0 $SEED & run Zn 1.0 0.0 $SEED & run Ca 0.0 1.0 $SEED & wait
  run MgZnCa 0.35 0.05 $SEED & 
done
wait
python3 postprocess_mixing.py > mixing_summary.txt
