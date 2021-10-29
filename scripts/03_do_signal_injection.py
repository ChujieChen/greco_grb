#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT: combo/stable
#!/usr/bin/env python
"""
1 GBM-GRB + 1 nonGBM-GRB = 2 GRBs (two different qsub)
6 time windows
40  different injections (0 ~ 9.1, step=0.2)
500 trials for each n_inj (run on 1 core)
20 s/trial
= 2k cpu hr = 1000 hivecpu * 2 hr

tws_in_second    = [  10,     25,    50,   100,   250,   500][::-1] 
"""

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
p = argparse.ArgumentParser(description="Signal Trials",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_name", default="GRB180423A", type=str, help="GRB name: GRByymmddC")
p.add_argument("--n_inj", default=0, type=float, help="Number of (poisson mean) injection")
p.add_argument("--n_trials", default=500, type=int, help="Number of trials")
p.add_argument("--tw_in_second", default=10, type=int, help="Length of the time window in seconds")
p.add_argument("--ncpu", default=4, type=int, help="Number of cores used")
p.add_argument("--batchIndex", default=0, type=int, help="Current batchIdx for this n_inj with this tw_in_second")
p.add_argument("--use_poisson", default=True, type=bool, help="Use poisson for n_inj")
args = p.parse_args()
###########################################################################

### testing on jupyter ###
# class args:
#     grb_name = "GRB180423A"    # real healpix example
#     # grb_name = "GRB190415A"    # fake healpix example
#     # grb_name = "GRB170529A"
#     n_inj = 10
#     n_trials = 10
#     tw_in_second = 10
#     ncpu = 4
#     batchIndex = 0
#     use_poisson = True
##########################

print("\n===== Loading no-healpix df =====\n")
df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")
df.head()

print("\n===== Loading healpix of {}=====\n".format(args.grb_name))   
try:
    healpix = np.load(DATA_DIR+"/grbwebgbm/healpix/{}_healpix_nside64.npy".format(args.grb_name))
    # healpix can contain negative values: faults due to Fermi-GBM
    healpix = np.maximum(healpix,0)
    ########## healpix reduce (< instead of <=) ##########
    healpix[healpix < isf_healpix(healpix, q=0.99)] = 0
except:
    raise Exception("Cannot load the healpix for grb: {}\n".format(args.grb_name))
    
grb_row = df.loc[df['grb_name'] == args.grb_name]
tw = args.tw_in_second/86400.
tw_start = grb_row.t_center - 0.5*tw
ra = grb_row.ra
dec = grb_row.dec
print("\n...Done\n")

print("\n===== Setting up csky =====\n")
data_dir = ICDATA_DIR
data_filenames = sorted(glob(data_dir + '/IC86_20*.data.npy'))
sig_filenames = sorted(glob(data_dir + '/IC86_2012.nu*_merged.npy'))
grl_filenames = sorted(glob(data_dir + '/GRL/IC86_20*.data.npy'))

################ energy lower bound #############
min_log_e = np.log10(10)
#################################################
bins_sindec = np.linspace(-1, 1, 25+1)  
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
                      , load_sig=True)  # false to save memory if needed 

#### used for spatial_prior_trial_runner
conf = {
    'ana': ana,
    #### llh basics: csky.conf
    'space': 'ps', # ps/fitps/template/prior
    'time': 'transient', # utf/lc/transient
    'energy': 'customflux', # fit/customflux
    'flux': cy.hyp.PowerLawFlux(2.5),
    #### inj.py - prior has some duplications against space's prior
    'sig': 'transient', # ps/tw/lc/transient/template/prior
    'full_sky': True,
    'extended': True,
    'mp_cpus': args.ncpu,
    'cut_n_sigma': 3
    }
cy.CONF.update(conf)

print("\n...Done\n")

print("\n===== Loading precomputed Bkg TS distribution =====\n")
bg_files = glob(ANA_DIR+"/allsky_scan/with_prior_background/tw{}/{}*.npz".format(args.tw_in_second, args.grb_name))
bg = cy.dists.Chi2TSD(np.ravel([sparse.load_npz(bg_file).toarray() for bg_file in bg_files]))
print("\n...Done\n")

print("\n===== Do trials =====\n")
src = cy.sources(
    ra=ra,
    dec=dec,
    deg=True,
    mjd=tw_start, 
    sigma_t=np.zeros_like(tw), 
    t_100=tw,  # in days
    prior=[hl.heal.HealHist(healpix)],
    name=args.grb_name
)

seed = abs(java_hash(src.name[0]+"_signal_batchIndex{}".format(args.batchIndex)))

"""
sptr = cy.get_spatial_prior_trial_runner(conf=cy.CONF
                                         ,src_tr=src
                                         ,llh_priors=[healpix]
                                         ,cut_n_sigma=5.) # src_tr is must for transient
"""

tr = cy.get_trial_runner(conf=cy.CONF, ana=ana, src=src)

with time('Doing injections'):
    """
    trials = sptr.get_many_fits(args.n_trials, 
                          n_sig=args.n_inj, 
                          poisson=args.use_poisson, 
                          seed=seed, 
                          logging=False)
    """
    trials = tr.get_many_fits(args.n_trials, 
                      n_sig=args.n_inj, 
                      poisson=args.use_poisson, 
                      seed=0, 
                      logging=False)
    
print("\n...Done\n")
print("\n===== Saving results =====\n")
outfilename = "{}_batchSize{}_batchIndex{}_tw{}_ninj{}.npy".format(args.grb_name, 
                                                                    args.n_trials, 
                                                                    args.batchIndex, 
                                                                    args.tw_in_second, 
                                                                   args.n_inj)

output_folder = cy.utils.ensure_dir(ANA_DIR+"/prior_injection/tw{}/{}".format(args.tw_in_second, 
                                                                              args.grb_name))
np.save(output_folder + "/" + outfilename, trials.as_array)
print("\nAll Done\n")







