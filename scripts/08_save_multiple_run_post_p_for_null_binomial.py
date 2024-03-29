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
p.add_argument("--batch_idx", default=0, type=int, help="starting index of runs")
args = p.parse_args()
###########################################################################

### testing on jupyter ###
# class args:
#     num_run = 3
#     batch_idx = 29
##########################

def get_post_p(pre_post, pre_p):
    """
    effective trial correction due to multiple time windows
    
    Parameters
    ----------
        pre_post: ndarray
            Npy from (ANA_DIR + f"/effective_trial/pre_post/{grb_name}_effective_trial.npy")
            dtype=[('pre_trial_p', '<f4'), ('post_trial_p', '<f4')])
        pre_p: float or array_like
            one or multiple pre_trial p value(s)
        
    Returns
    -------
        post_p: float or array_like
            one or multiple post_trial p value(s)
    """
    pre_p = np.float32(pre_p)
    idx = np.searchsorted(pre_post['pre_trial_p'], pre_p, side='right') - 1
    idx = np.maximum(idx, 0)
    idx = np.minimum(idx, pre_post['post_trial_p'].shape[0] - 1)
    return pre_post['post_trial_p'][idx]

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

def save_multiple_run_post_p(num_run, batch_idx=0):
    multiple_run_post_p = get_multiple_run_post_p(num_run, batch_idx=batch_idx, load=False)
    np.save(ANA_DIR+f"/binomial_test/null_binom/multiple_run_post_p/multiple_run_post_p_numRun{num_run}_batchIdx{batch_idx}.npy", 
            multiple_run_post_p)
    
with time("test get_multiple_run_post_p"):
    save_multiple_run_post_p(args.num_run, args.batch_idx)
    
    