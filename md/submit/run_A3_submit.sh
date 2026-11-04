#!/bin/bash
# A3 submit helper -- one 1-core job per account, from the folder that holds the case folders and detect_shear_bands.py.
#   p230674mt (W2_new/new_md_package):   bash run_A3_submit.sh glass0_r5 glass0_r1 glass0_r6 vf_r5_1 vf_r5_2 vf_r5_3 temp_350 temp_400 rate2e9_run4
#   p240165mt (.../codes/new_md_package): bash run_A3_submit.sh rep_r5_a rep_r5_b rep_r5_c rep_r1_a rep_r1_b rep_r1_c rep_r6_a rep_r6_b rep_r6_c temp_250 temp_300 rate2e9_run3
# Uses 1 queue slot.  Partition/time can be overridden: PART=shortq TIME=04:00:00 bash run_A3_submit.sh ...
PART=${PART:-mediumq}; TIME=${TIME:-04:00:00}
sbatch --partition=$PART --time=$TIME --ntasks=1 --cpus-per-task=4 --job-name=A3_shear --output=A3_%j.out --error=A3_%j.err \
  --wrap="source ~/shap_env/bin/activate; cd $PWD; python3 run_A3_shearbands.py $*"
