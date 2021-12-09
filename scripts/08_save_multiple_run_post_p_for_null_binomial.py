import os
import sys
import numpy as np
import healpy as hp
import histlite as hl
import csky as cy
import pandas as pd
from scipy import sparse

from glob import glob
timer = cy.timing.Timer()
time = timer.time

sys.path.append('../../')
from greco_grb.scripts import SETTING
paths = SETTING.PATH()
print(paths)
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR

from greco_grb.scripts.utils import *

import argparse
######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Save multiple background p-values(effective corrected) for null hypothesis binomial test",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--num_run", default=1, type=int, help="number of runs")
args = p.parse_args()
###########################################################################

def get_post_p(grb_name, pre_p):
    """
    effective trial correction due to multiple time windows
    
    Parameters
    ----------
        grb_name: str
            name of the grb
        pre_p: float or array_like
            one or multiple pre_trial p value(s)
        
    Returns
    -------
        post_p: float or array_like
            one or multiple post_trial p value(s)
    """
    pre_p = np.float32(pre_p)
    pre_post = np.load(ANA_DIR + f"/effective_trial/pre_post/{grb_name}_effective_trial.npy")
    idx = np.searchsorted(pre_post['pre_trial_p'], pre_p, side='right') - 1
    idx = np.maximum(idx, 0)
    idx = np.minimum(idx, pre_post['post_trial_p'].shape[0] - 1)
    return pre_post['post_trial_p'][idx]

import multiprocessing as mp

def get_grb_pre_post_p(grb_name, trial_no):
    """
    A helper function for function:get_all_GRB_best_p_values
    """
    pre_trial_tw_p = np.load(ANA_DIR + f"/effective_trial/min_tw_p/{grb_name}_min_tw_p.npy")[trial_no]
    post_trial_p = get_post_p(grb_name, pre_trial_tw_p[1])
    return (grb_name, pre_trial_tw_p[0], pre_trial_tw_p[1], post_trial_p)

def get_all_GRB_best_p_values(trial_no):
    """
    Get 2268 p-values at the trial_no-th trial. 
    An array of (grb_name, tw_index, pre_p, post_p) will be returned
    The post_p are p-values corrected with `effective_trial/pre_post`
    
    Parameters
    ----------
        trial_no: int
            Some value between 0 and N (N is 1,000,000)
        
    Returns
    -------
        grb_tw_pre_post_p: array_like
            2268 (grb_name, tw_index, pre_p, post_p)
    """
        
    df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")
    grb_tw_pre_post_p = []
    with mp.Pool() as pool:
        grb_tw_pre_post_p = np.array(
            pool.starmap(
                get_grb_pre_post_p, 
                zip(list(df.grb_name.values),[trial_no]*len(df.grb_name))
            )
        )
    dt = np.dtype([('grb_name', 'U10'), 
                   ('tw_idx', np.intc), 
                   ('pre_trial_p', np.float32), 
                   ('post_trial_p', np.float32)])
    grb_tw_pre_post_p = np.array(
        list(map(tuple, grb_tw_pre_post_p)), 
        dtype=dt
    )
    return grb_tw_pre_post_p
    
# average null hypothesis Binomial

def get_multiple_run_post_p(num_run):
    """
    Get num_run * 2268 p-values. 
    The post_p are p-values corrected with `effective_trial/pre_post`
    
    Parameters
    ----------
        num_run: int
            First a few runs (e.g. 1000)
        
    Returns
    -------
        multiple_run_post_p: array_like
            shape: num_run * 2268
    """
    multiple_run_post_p = []
    # get_all_GRB_best_p_values is using multiprocessing
    # so we cannot use multiprocessing here
    for i in range(num_run):
        multiple_run_post_p.append(get_all_GRB_best_p_values(i)['post_trial_p'])
    multiple_run_post_p = np.array(multiple_run_post_p)
    return multiple_run_post_p

def save_multiple_run_post_p(num_run):
    multiple_run_post_p = get_multiple_run_post_p(num_run)
    np.save(ANA_DIR+f"/binomial_test/null_binom/multiple_run_post_p_{num_run}.npy", 
            multiple_run_post_p)
with time(f"Run {args.num_run} runs"):   
    save_multiple_run_post_p(args.num_run)