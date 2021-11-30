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
p = argparse.ArgumentParser(description="Save effective trial mapping: pre-trial p vs. post-trial p",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_idx_start", default=0, type=int, help="GRB idx start: [start, end)")
p.add_argument("--grb_idx_end", default=1, type=int, help="GRB idx end: [start, end)")
args = p.parse_args()
###########################################################################


def save_effective_trial(grb_name):
    bg_files = sorted(glob(ANA_DIR + f"/allsky_scan/with_prior_background/tw*/{grb_name}_*.npz"), 
                 key=lambda x: int(x[x.find("/tw")+3:x.find(f"/{grb_name}")]))
    assert len(bg_files) > 0, print(f"cannot find bkg files of {grb_name}")

    bgs = np.array([sparse.load_npz(bg_file).toarray()[0] for bg_file in bg_files])
    bgs_sorted = np.apply_along_axis(sorted, 1, bgs)

    pvals = []
    for bg, bg_sorted in zip(bgs, bgs_sorted):
        pvals.append(np.apply_along_axis(lambda x: (bg.size - np.searchsorted(bg_sorted, x, side='left')) / bg.size, 0, bg))
    pvals = np.array(pvals)
    best_pvals = np.sort(pvals.min(axis=0))

    hist, bin_edges = np.histogram(best_pvals, 
                               bins=np.r_[np.unique(best_pvals),1.2], 
                               density=True)
    pre_trial_p = bin_edges
    post_trial_p = np.r_[np.cumsum(hist*np.diff(bin_edges)), 1]
    dt = np.dtype([('pre_trial_p', np.float32), ('post_trial_p', np.float32)])

    pre_post = np.transpose(np.array([pre_trial_p, post_trial_p]))
    pre_post = np.array(list(map(tuple, pre_post)), dtype=dt)

    np.save(ANA_DIR + f"/effective_trial/{grb_name}_effective_trial.npy", pre_post)
    
    
df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")
for i in range(args.grb_idx_start, args.grb_idx_end):
    try:
        grb_name = df.grb_name.iloc[i]
        save_effective_trial(grb_name)
    except:
        raise Exception("Something wrong with grb idx {}: {}\n".format(i, grb_name))