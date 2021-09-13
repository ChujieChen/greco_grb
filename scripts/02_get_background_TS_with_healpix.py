import os
import sys
print("Python version: ", end=' ')
print(sys.version)

import numpy as np
import healpy as hp
# import histlite as hl
import csky as cy
import pandas as pd
from scipy import sparse

# import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as colors
############# comment out below two lines on clusters ##########
# %matplotlib inline
# %matplotlib notebook
################################################################
from glob import glob
timer = cy.timing.Timer()
time = timer.time

###### Local Import ######
import SETTING
paths = SETTING.PATH(osg=False)
print(paths)
LOCATION = paths.LOCATION
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR

from utils import *

import argparse
######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Background Trials",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_name_idx_start", default=0, type=int, help="Starting index of one GRB")
p.add_argument("--grb_name_idx_end", default=1, type=int, help="Ending index of one GRB")
p.add_argument("--tw_in_second", default=10, type=int, help="Length of the time window in seconds")
args = p.parse_args()
###########################################################################


### testing on jupyter ###
# class args:
#     grb_name_idx_start = 0    # [start, end)   [inclusive, exclusive)
#     grb_name_idx_end = 2
#     # grb_name = "GRB180423A"    # real healpix example
#     # grb_name = "GRB190611B"    # fake healpix example
#     tw_in_second = 25
##########################

print("\n===== Loading no-healpix df =====\n")
df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHeaplix_2297.pkl")
print("\n===== Done =====\n")

print("\n===== Loading no-healpix background all-sky scan =====\n")
files = glob(ANA_DIR + "/allsky_scan/no_prior_versatile/tw{}/GRB180423A_batchSize*_batchIndex*_tw{}.npz".format(args.tw_in_second, args.tw_in_second))
files.sort(key=lambda x: int(x[x.find("batchIndex") + len("batchIndex"): x.find("_tw")]))
print("\n===== Done =====\n")

def getCombinedMaxTS(healpix, scan):
    """
    Combine the non-prior background TS scan with the healpix prior
    then return the maximum TS
    
    Parameters
    ----------
    healpix : array_like
        A prior healpix map. shape [hp.nside2pix(64), ]
        
    scan : array_like (scipy.sparse matrix)
        One non-prior background scan (shape [1, hp.nside2pix(64)])
    
    Returns
    -------
    TS : scalar
        The maximum TS from the combined background.
    
    See also
    -------
    N/A
    """
    healpix = np.maximum(healpix, 1e-15)
    ts = scan.copy()
    # non-overlapping region will have negative TS anyways
    ts.data += 2. * (np.log(healpix[scan.indices]) - np.log(np.max(healpix)))
    ts.data[~np.isfinite(ts.data)] = 0
    return np.max(ts)


for grb_idx in range(args.grb_name_idx_start, args. grb_name_idx_end):
    grb_name = df.grb_name[grb_idx]
    print("\n===== Loading Prior {} =====\n".format(grb_name))
    probs = np.maximum(1e-15, np.load(DATA_DIR+"/grbwebgbm/healpix/{}_healpix_nside64.npy".format(grb_name)))
    print("\n===== Done =====\n")
    
    print("\n===== Calculating MAXs =====\n")
    with time("{}_tw{}".format(grb_name, args.tw_in_second)):
        TSs = []
        for i, file in enumerate(files):
            if i % (len(files)/10) == 0:
                print("Working on number of file: {}".format(i))
            scans = sparse.load_npz(file)
            for idx in range(scans.shape[0]):
                scan = scans[idx]
                TS = getCombinedMaxTS(probs, scan)
                TSs.append(TS)
    print("\n===== Done =====\n")
    
    
    print("\n===== Saving... =====\n")
    sTSs = sparse.csr_matrix(TSs, dtype=float)
    
    outfilename = "{}_tw{}_NTrial{}.npz".format(grb_name, 
                                                 args.tw_in_second,
                                                sTSs.shape[1])
    output_folder = cy.utils.ensure_dir(ANA_DIR+"/allsky_scan/with_prior_background/tw{}".format(args.tw_in_second))
    sparse.save_npz(output_folder+"/{}".format(outfilename) ,sTSs)
    print("\n===== Saved =====\n")
    
print("All done")