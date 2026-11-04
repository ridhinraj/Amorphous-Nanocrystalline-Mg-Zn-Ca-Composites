#!/bin/bash
#SBATCH -p mediumq
#SBATCH -n 16
#SBATCH -t 06:00:00
#SBATCH -J mixfix
cd /gpfs-home/p230674mt/W2_new/new_md_package/checks
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
LMP=/gpfs-home/p230674mt/lammps_src/src/lmp_mpi
run() { mpirun -np 4 $LMP -in in.mixing_check -var TAG $1 -var XZN $2 -var XCA $3 -var SEED $4 -var LAT $5 -var TMELT $6 -var ONLY1000 $7 -screen none; }
# Zn at 1000 K only (new protocol), seeds 1-3, plus the Mg run that was killed
run Zn 1.0 0.0 1 3.95 1000.0 1 & run Zn 1.0 0.0 2 3.95 1000.0 1 & run Zn 1.0 0.0 3 3.95 1000.0 1 & run Mg 0.0 0.0 1 4.5 2000.0 0 & wait
# the other runs that were killed by the out-of-memory event (skip any whose log already ends in DONE)
run Mg 0.0 0.0 3 4.5 2000.0 0 & run Ca 0.0 1.0 3 4.5 2000.0 0 & run MgZnCa 0.35 0.05 2 4.5 2000.0 0 & wait
python3 postprocess_mixing.py > mixing_summary.txt
