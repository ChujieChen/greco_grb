#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT: combo/stable
#!/usr/bin/env python
r"""
Perform all-sky scans with no spatial prior information.
Record the TS (and number of true and fitted events around the
best fit location) for each pixel, creating a background-only TSD
at each point in the sky.
Combine this map with GBM healpix maps or grab a single pixel
for non-GBM bursts.
"""
###################### OSG Specific ###########################
# import sys
# modules_dir = '/cvmfs/icecube.opensciencegrid.org/users/cjchen'
# sys.path.append(modules_dir+'/csky')
# sys.path.append(modules_dir+'/greco_grb')
# sys.path.append(modules_dir+'/python3.7/site-packages/')
###############################################################
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

tw_batchSize_map = {
                 10:1000,       # 59 cpu seconds/trial
                 25:1000,       # 55 cpu seconds/trial
                 50:1000,       # 54 cpu seconds/trial
                 100:1000,      # 56 cpu seconds/trial
                 250:1000,      # 66 cpu seconds/trial
                 500:1000       # 67 cpu seconds/trial
                }

import argparse

######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Background Trials",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_name", default="GRB180423A", type=str, help="Name of one GRB")
p.add_argument("--batchNtrials", default=10, type=int, help="Number of trials in this batch")
p.add_argument("--batchIndex", default=0, type=int, help="Index of current batch")
p.add_argument("--tw_in_second", default=10, type=int, help="Length of the time window in seconds")
p.add_argument("--ncpu", default=1, type=int, help="Number of CPU to give Csky")
p.add_argument("--mode", default="production", type=str, help="Mode: production or testing")
p.add_argument("--outfilename", default="", type=str, help="Output filename should have type .npz. Highly recommended on OSG")
args = p.parse_args()
###########################################################################

### testing on jupyter ###
"""
class args:
    grb_name = "GRB180423A"    # real healpix example
    # grb_name = "GRB190415A/GRB190611B"    # fake healpix example
    batchNtrials = 40
    batchIndex = 0
    tw_in_second = 10
    ncpu = 4
    mode = "testing"
    outfilename = ""
"""
##########################
if tw_batchSize_map[args.tw_in_second] != args.batchNtrials and args.mode=="prodection":
    raise Exception("You are not using the recommended batchNtrials wrt this tw_in_second!")
    
print("\n===== Loading noHeaGRB list =====\n") 
if args.mode != "production" and LOCATION != "IC-OSG":
    ANA_DIR = cy.utils.ensure_dir(ANA_DIR + "/test")
# All times in days, all angles in degrees
try:
    print("Loading pkl in DATA_DID")
#     df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHeaplix.pkl")
    df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHeaplix_2297.pkl")
    
except:
    try:
        print("Loading pkl in current path")
#         df = pd.read_pickle("grbweb_gbm_noHeaplix.pkl")
        df = pd.read_pickle("grbweb_gbm_noHeaplix_2297.pkl")
    
    except:
        raise Exception("Cannot pd.reade_picle() the grbweb_gbm_noHeaplix.pkl.\n")

print("\n...Done\n")
    
grb_row = df.loc[df['grb_name'] == args.grb_name]
tw = args.tw_in_second/86400.
tw_start = grb_row.t_center - 0.5*tw
ra = grb_row.ra
dec = grb_row.dec
print("\n...Done\n")

print("\n===== Setting up csky =====\n")
data_dir = ICDATA_DIR
# data_filenames = sorted(glob(data_dir + '/IC86_20*.data_with_angErr.npy'))
data_filenames = sorted(glob(data_dir + '/IC86_20*.data.npy'))
# sig_filenames = sorted(glob(data_dir + '/IC86_2012.nu*_merged_with_angErr.npy'))
# load nue only to save memory, never used in this .py
sig_filenames = sorted(glob(data_dir + '/IC86_2012.nu*_merged.npy'))
grl_filenames = sorted(glob(data_dir + '/GRL/IC86_20*.data.npy'))

################ energy lower bound #############
min_log_e = np.log10(10)
#################################################
bins_sindec = np.linspace(-1, 1, 25+1)  
# bins_logenergy = np.linspace(min_log_e, 5, 30+1)
bins_logenergy = np.linspace(min_log_e, 4, 25+1)

data = [np.load(data_filename) for data_filename in data_filenames]
data = np.hstack(data)
sig = [np.load(sig_filename) for sig_filename in sig_filenames]
sig = np.hstack(sig)
grl = [np.load(grl_filename) for grl_filename in grl_filenames]
grl = np.hstack(grl)
if min_log_e is not None:
    data_mask = data['logE'] > min_log_e
    data = data[data_mask]
    sig_mask = sig['logE'] > min_log_e
    sig = sig[sig_mask]
    
dataset_spec = cy.selections.CustomDataSpecs.CustomDataSpec(data, sig, np.sum(grl['livetime']),
                                                     sindec_bins=bins_sindec,
                                                     logenergy_bins=bins_logenergy,
                                                     grl=grl, key='greco_v2.10', cascades=True)

ANA_DIR = cy.utils.ensure_dir(ANA_DIR)
# on OSG
# ana_dir = "./"
ana = cy.get_analysis(cy.selections.repo
                      , dataset_spec
                      , dir=ANA_DIR
                      , load_sig=False)  # to save memory    

conf = {
    'ana': ana,
    #### llh basics: csky.conf
    'space': 'ps', # ps/fitps/template/prior
    'time': 'transient', # utf/lc/transient
    'energy': 'customflux', # fit/customflux
    'flux': cy.hyp.PowerLawFlux(2),
    #### inj.py - prior has some duplications against space's prior
    'sig': 'transient', # ps/tw/lc/transient/template/prior
    'full_sky': True,
    'extended': True,
    'mp_cpus': args.ncpu,
    'cut_n_sigma': 3
    }

cy.CONF.update(conf)

print("\n===== Generating seeds for current batch =====\n")   
src = cy.utils.Sources(
    ra=ra,
    dec=dec,
    deg=True,
    mjd=tw_start, 
    sigma_t=np.zeros_like(tw), 
    t_100=tw,  # in days
    # prior=[hl.heal.HealHist(healpix)],
    name=args.grb_name
)
sstr = cy.get_sky_scan_trial_runner(conf=cy.CONF
                                    ,nside=64
                                    ,src_tr=src)

rng=np.random.default_rng(abs(java_hash(src.name[0])))
seeds = rng.integers(int(1e9), size=int(2e8))[args.batchNtrials*args.batchIndex: args.batchNtrials*(args.batchIndex + 1)]
print("\n...Done\n")

print("\n===== Scanning =====\n") 
all_pixel_TS = sparse.lil_matrix((len(seeds), hp.nside2npix(64)), dtype=float)
with time("allsky scramble scan"):
    for no_trial, seed in enumerate(seeds):
        if no_trial % (len(seeds) // 10) == 0:
            print("Working on no_trial: {} \n".format(no_trial))
        # scan (3,49152): -log10p, TS, ns
        scan = sstr.get_one_scan(n_sig=0
                                 , poisson=False
                                 , seed=seed
                                 , TRUTH=False
                                 , mp_cpus=args.ncpu
                                 , logging=False)
        all_pixel_TS[no_trial] = scan[1]
print("\n...Done\n")

print("\n===== Converting to scipy.sparse and Save to disk =====\n")
with time("To scipy.sparse npz"):
    hp_sparse = all_pixel_TS.tocsr()
    outfilename = args.outfilename
    if not outfilename:
        outfilename = "{}_batchSize{}_batchIndex{}_tw{}.npz".format(args.grb_name, 
                                                                    args.batchNtrials, 
                                                                    args.batchIndex, 
                                                                    args.tw_in_second)
## on locations other than OSG
    output_folder = cy.utils.ensure_dir(ANA_DIR+"/allsky_scan/no_prior_versatile/tw{}".format(args.tw_in_second))
    sparse.save_npz(output_folder+"/{}".format(outfilename)
                    ,hp_sparse)
## on OSG
#     sparse.save_npz("{}".format(outfilename)
#                     ,hp_sparse)           
##
print("######## All Done. ###########")