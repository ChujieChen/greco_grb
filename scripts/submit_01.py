import sys
import os
try:
    from os import errno
except:
    import errno
import glob
import numpy as np
import pandas as pd

####### Check before running 1/3 ########
N_TRIALS = int(2e3)
MODE = "testing" # production or testing
#######################################

def ensure_dir(dirname):
    """Make sure ``dirname`` exists and is a directory."""
    if not os.path.isdir(dirname):
        try:
            os.makedirs(dirname)   # throws if exists as file
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    return dirname

print("\n===== Making scratch folders =====\n")
jobdir =  '/storage/home/hhive1/cchen641/scratch/icecube/job_pbs_out'
ensure_dir(jobdir)
pbs_dir = jobdir+"/pbs"
ensure_dir(pbs_dir)
out_dir = jobdir+"/out"
ensure_dir(out_dir)
print("\n=====  Done  =====\n")

print("\n===== Getting Paths =====\n")
import SETTING
paths = SETTING.PATH()
print(paths)
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR
print("\n=====  Done  =====\n")

print("\n===== Loading GRB list =====\n")
try:
    df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHeaplix.pkl")
except:
    raise Exception("Cannot pd.reade_pickle() the grbweb_gbm_noHeaplix.pkl.\n")
print("\n=====  Done  =====\n")

print("\n===== Creating PBS Scripts =====\n")
########### Check before running 2/3 ##############
tw_batchSize_map = {1:1000000,
                2:1000000,
                5: 1000,  # 5:1000000,
                 10:1000000,
                 20:100000,
                 50: 1000, # 50:100000,
                 100:100000,
                 200:100000,
                 500: 1000, # 500:100000,
                 1000:100000,
                 2000:100000,
                 5000: 1000, # 5000:100000 
                }
####################################################
########### Check before running 3/3 ##############
# tws_in_second = [1,2,5,10,20,50,100,200,500,1000,2000,5000]
tws_in_second = [5,50,500,5000]  # for testing only
# grbs = df.grb_name
# grbs = df.grb_name[2048:2049]  # for testing only: real healpix: GRB120427B
grbs = df.grb_name[3:4]        # for testing only: fake healpix: GRB190611B
#######################################################
for tw_in_second in tws_in_second:
    batchSize = int(tw_batchSize_map[tw_in_second])
    print("\n===== TW = {} second, batchSize = {} trials =====\n".format(tw_in_second, batchSize))
    
    for grb in grbs:
        print("\n===== {} =====\n".format(grb))
        pbs_tw_grb_dir = ensure_dir(pbs_dir+"/{}/tw{}".format(grb, tw_in_second))
        
        for idx in range(N_TRIALS // batchSize):
            if idx % (N_TRIALS // batchSize / 10 if N_TRIALS // batchSize > 10 else 1) == 0:
                print("\n===== batchIndex = {} =====\n".format(idx))
            pbsname = pbs_tw_grb_dir + "/{}_tw{}_batchSize{}_batchIndex{}.pbs".format(grb, tw_in_second, batchSize, idx)
            script = open(pbsname, "w") 
            out_tw_grb_dir = ensure_dir(out_dir+"/{}/tw{}".format(grb, tw_in_second))
            outname = out_tw_grb_dir + "/{}_tw{}_batchSize{}_batchIndex{}.out".format(grb, tw_in_second, batchSize, idx)
            pythonExec = "python 01_do_background_trials.py --grb_name {} --tw_in_second {} --batchNtrials {} --batchIndex {} --mode {}".format(grb, tw_in_second, batchSize, idx, MODE)
            if idx == (N_TRIALS // batchSize - 1):
                pythonExec = "python 01_do_background_trials.py --grb_name {} --tw_in_second {} --batchNtrials {} --batchIndex {}  --concat 1 --totalNtrials {} --mode {}".format(grb, tw_in_second, batchSize, idx, N_TRIALS, MODE)
                
            delete_out = "rm {}".format(outname)   
            if MODE == "testing":
                delete_out = " "
            exeScript = [
                "#PBS -N {}_tw{}_batchSize{}_batchIndex{}".format(grb, tw_in_second, batchSize, idx),
                "#PBS -l nodes=1:ppn=1",
                "#PBS -l pmem=1gb",
                "#PBS -l walltime=3:00:00",
                "#PBS -q hive",
                "#PBS -j oe",
                "#PBS -o {}".format(outname),
                " ",
                "cd /storage/home/hhive1/cchen641/icecube/greco_grb/scripts",
                "icpy3",
                pythonExec,
                delete_out  # /scratch has a limit of 1 million files
            ]
            # No worries: this will erase the previous content and write content
            for line in exeScript:
                script.write("%s\n" % line) 
            script.close()
            
            cmd = "qsub -q hive {}".format(pbsname)
            os.system(cmd)
            
            cmd = "rm {}".format(pbsname) # /scratch has a limit of 1 million files
            os.system(cmd)
            
