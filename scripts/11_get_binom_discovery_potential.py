#!/usr/bin/env python
import os
import sys
print("Python version: ", end=' ')
print(sys.version)

import numpy as np
import healpy as hp
import histlite as hl
import csky as cy
import pandas as pd
from scipy import sparse

############# comment out below lines on clusters ##########
# import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as colors
# %matplotlib inline
# %matplotlib notebook
############################################################

from glob import glob
timer = cy.timing.Timer()
time = timer.time

###### Local Import ######
sys.path.append('../../')
from greco_grb.scripts import SETTING
paths = SETTING.PATH()
print(paths)
LOCATION = paths.LOCATION
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR
from greco_grb.scripts.utils import *

import argparse
######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Get multiple background binomial test results",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--seed", default=0, type=int, help="seed for the permutation")
args = p.parse_args()
###########################################################################
### testing on jupyter ###
# class args:
#     seed = 0
##########################



from scipy import stats
def binomial_test(p_values, epsilon=None):
    """
    Perform IceCube binomial test
    
    Parameters
    ----------
        p_values: array_like
            an array of p-values
        
    Returns
    -------
        best_k: int
            number of p-values that minimizes the binomial probability
            Note this is 1-index'ed
        
        best_binomial_p: float
            the optimized binomial probability
            
        threshold_p_value: float
            the best_k-th p_value (a.k.a. p_k)
            this GRB and all GRBs having p_values smaller than this threshold_p_value are important
            
        indices: array_like
            GRB indices after sorting (assume the input p_values' indices are from 1 - 2268)
            
        
        binomial_ps: array_like
            binomial probabilities for different k=1,2,...,N
    """
    p_values = np.array(p_values)
    N = p_values.size
    ks = np.argsort(p_values)
    p_values = np.sort(p_values)
    binomial_ps = stats.binom.sf(np.r_[0:N], N, p_values)
    best_k = np.argmin(binomial_ps) + 1
    best_binomial_p = binomial_ps[best_k - 1]
    if (epsilon is not None) and (1.0 - best_binomial_p < epsilon):
        return N, 1, p_values[-1], ks, binomial_ps
    threshold_p_value = p_values[best_k - 1]
    return best_k, best_binomial_p, threshold_p_value, ks, binomial_ps

def get_multiple_run_post_p(num_run, batch_idx=0, load=False):
    """
    Get num_run * 2268 p-values. 
    The post_p are p-values corrected with `effective_trial/pre_post`
    
    Parameters
    ----------
        num_run: int
            First a few runs starting from batch_idx * num_run (e.g. 1000)
        
        batch_idx: int
            get num_run runs from start_idx (e.g.  [batch_idx*num_run, (batch_idx+1)*num_run))
            if -1 (only when load=True), get all batches with num_run
            
        load: bool
            are we loading or making runs
        
    Returns
    -------
        multiple_run_post_p: array_like
            shape: num_run * 2268
    """
    print(f"Running num_run {num_run}, batch_idx {batch_idx}")
    multiple_run_post_p = []
    # get_all_GRB_best_p_values is using multiprocessing
    # so we cannot use multiprocessing here
    if batch_idx==-1:
        assert load == True
    if load:
        if batch_idx==-1:
            files = glob(ANA_DIR+f"/binomial_test/null_binom/multiple_run_post_p/multiple_run_post_p_numRun{num_run}_batchIdx*.npy")
            files = sorted(files, key=lambda x: int(x[x.find("_batchIdx")+9:x.find(".npy")]))
        else:
            files= glob(ANA_DIR+f"/binomial_test/null_binom/multiple_run_post_p/multiple_run_post_p_numRun{num_run}_batchIdx{batch_idx}.npy")
        return np.vstack([np.load(x) for x in files])
    
    # load all GRB names
    df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")
    for grb_name in df.grb_name.values:
        # load pre_trial_tw_p_grb for this grb
        pre_trial_tw_p_grb = np.load(ANA_DIR + f"/effective_trial/min_tw_p/{grb_name}_min_tw_p.npy")
        # load for pre_post this grb
        pre_post = np.load(ANA_DIR + f"/effective_trial/pre_post/{grb_name}_effective_trial.npy")
        
        ## multiple_run_post_p_grb for this grb
        multiple_run_post_p_grb = []
        for i in range(batch_idx*num_run, (batch_idx+1)*num_run):
            pre_trial_tw_p = pre_trial_tw_p_grb[i]
            single_run_post_p_grb = get_post_p(pre_post, pre_trial_tw_p[1])
            multiple_run_post_p_grb.append(single_run_post_p_grb)
        multiple_run_post_p.append(multiple_run_post_p_grb)
        
    # make the return to shape N*2268
    multiple_run_post_p = np.array(multiple_run_post_p).transpose()
    return multiple_run_post_p



"""
Note: grb_name indices will NOT get permuted
"""
## all 1,000,000 runs, load=True ONLY
multiple_run_post_p = get_multiple_run_post_p(2500, batch_idx=-1, load=True)

seed = args.seed
print(f"seed is {seed}")

rng = np.random.default_rng(seed)
with time(f"Permute post p's"):
    multiple_run_post_p_permuted = rng.permuted(multiple_run_post_p, axis=0)

print(f"Running binomial test with seed {seed}")
with time(f"Get binomial test for {multiple_run_post_p.shape[0]} runs"):
    multiple_permuted_null_binomial_results = np.array([binomial_test(single_run_post_p, epsilon=None)[:3] for single_run_post_p in multiple_run_post_p_permuted])

print("Saving...")    
try:
    np.save(
        ANA_DIR+f"/binomial_test/null_binom/binomial_results/null_binomial_results_numRun{multiple_run_post_p.shape[0]}_seed{seed}.npy", 
        multiple_permuted_null_binomial_results
    )
except:
    print("An exception occurred when saving the result")

print("DONE")
