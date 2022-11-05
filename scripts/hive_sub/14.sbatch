#!/bin/bash#SBATCH -Jjobname=14
#SBATCH --account=hive-itaboada3
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=110:00:00
#SBATCH -phive
#SBATCH -oReport-%j.out
#SBATCH --array=1-378%400
cd /storage/home/hhive1/cchen641/icecube/greco_grb/scripts
icpy3
params=`sed -n "${PBS_ARRAYID} p" hive_sub/14.txt`
paramsArray=($params)
arg0=${paramsArray[0]}
arg1=${paramsArray[1]}
arg2=${paramsArray[2]}
arg3=${paramsArray[3]}
arg4=${paramsArray[4]}
srun python 14_stack_TS_background.py --grb_idx_start $arg0 --grb_idx_end $arg1 --batchNtrials $arg2 --batchIndex $arg3 --ncpu $arg4 