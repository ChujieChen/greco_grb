import os
import sys
import numpy as np
import healpy as hp
import histlite as hl
import csky as cy
import pandas as pd

# import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as colors
# %matplotlib inline
# %matplotlib notebook

from glob import glob
timer = cy.timing.Timer()
time = timer.time

import SETTING
paths = SETTING.PATH()
print(paths)

USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR

tw_batchSize_map = {1:1000000,
                2:1000000,
                5:1000000,
                 10:1000000,
                 20:100000,
                 50:100000,
                 100:100000,
                 200:100000,
                 500:100000,
                 1000:100000,
                 2000:100000,
                 5000:100000
                }

import argparse

######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Background Trials",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_name", default="GRB190612A", type=str, help="Name of one GRB")
p.add_argument("--batchNtrials", default=100, type=int, help="Number of trials in this batch")
p.add_argument("--batchIndex", default=0, type=int, help="Index of current batch")
p.add_argument("--tw_in_second", default=10, type=int, help="Length of the time window in seconds")
p.add_argument("--concat", default=0, type=int, help="True(1) for the last batchIndex. used to clear up everything")
p.add_argument("--totalNtrials", default=100000000, type=int, help="Number of total trials (1e8)")
p.add_argument("--ncpu", default=1, type=int, help="Number of CPU to give Csky")
p.add_argument("--mode", default="production", type=str, help="Mode: production or testing")
args = p.parse_args()
###########################################################################
if tw_batchSize_map[args.tw_in_second] != args.batchNtrials and args.mode=="prodection":
    raise Exception("You are not using the recommended batchNtrials wrt this tw_in_second!")

print("\n===== Loading GRB list =====\n") 
if args.mode != "production":
    ANA_DIR = ANA_DIR + "/test"
# All times in days, all angles in degrees
try:
    df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHeaplix.pkl")
except:
    raise Exception("Cannot pd.reade_picle() the grbweb_gbm_noHeaplix.pkl.\n")

print("\n...Done\n")
print("\n===== Loading healpix of {}=====\n".format(args.grb_name))   
try:
    healpix = np.load(DATA_DIR+"/grbwebgbm/healpix/{}_healpix_nside64.npy".format(args.grb_name))
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
data_filenames = sorted(glob(data_dir + '/IC86_20*.data_with_angErr.npy'))
sig_filenames = sorted(glob(data_dir + '/IC86_2012.nu*_merged_with_angErr.npy'))
grl_filenames = sorted(glob(data_dir + '/GRL/IC86_20*.data.npy'))

min_log_e = 1.0
bins_sindec = np.linspace(-1, 1, 25+1)  
bins_logenergy = np.linspace(min_log_e, 5, 30+1)

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
                                                     grl=grl, key='greco_v2.4', cascades=True)

ana_dir = cy.utils.ensure_dir(ANA_DIR)
ana = cy.get_analysis(cy.selections.repo, dataset_spec, dir=ana_dir, load_sig=True)    
conf = {
    'ana': ana,
    #### llh basics: csky.conf
    'space': 'prior', # ps/fitps/template/prior
    'time': 'transient', # utf/lc/transient
    'energy': 'customflux', # fit/customflux
    'flux': cy.hyp.PowerLawFlux(2),
    #### inj.py - prior has some duplications against space's prior
    'sig': 'transient', # ps/tw/lc/transient/template/prior
    'extended': True,
    }

cy.CONF['mp_cpus'] = args.ncpu
cy.CONF.update(conf)
    
print("\n===== Generating seeds for current batch =====\n")   
src = cy.utils.Sources(
    ra=ra,
    dec=dec,
    deg=True,
    mjd=tw_start, 
    sigma_t=np.zeros_like(tw), 
    t_100=tw,  # in days
    prior=[hl.heal.HealHist(healpix)],
    name=args.grb_name
)
tr = cy.get_trial_runner(src=src)

rng=np.random.default_rng(abs(hash(src.name[0])))
seeds = rng.integers(1e9, size=int(2e8))[args.batchNtrials*args.batchIndex: args.batchNtrials*(args.batchIndex + 1)]
print("\n...Done\n")
print("\n===== Getting fits =====\n")   
fits = []
with time('fits from n_trials'):
    for no_trial, seed in enumerate(seeds):
        if no_trial % (args.batchNtrials // 10) == 0:
            print("Working on no_trial: {} \n".format(no_trial))
        bg_trial = tr.get_one_trial(0, seed=seed)
        old_stdout = sys.stdout # suppress output start
        sys.stdout = open(os.devnull, "w")
        fit = tr.get_one_fit_from_trial(
            bg_trial, 
            TRUTH=False, 
            logging=False
        )
        sys.stdout = old_stdout # suppress output end
        fits.append(fit)
print("\n...Done\n")         
print("\n===== Saving fits =====\n")         
fits_df = pd.DataFrame(fits, columns=['TS', 'ns'])
sfits_df = fits_df.astype(pd.SparseDtype("float", 0.0))
output_folder = cy.utils.ensure_dir(ANA_DIR+"/bg_trials/{}/tw{}".format(args.grb_name, args.tw_in_second))
sfits_df.to_pickle(output_folder+"/{}_batchSize{}_batchIndex{}_tw{}.pkl".format(
    args.grb_name, 
    args.batchNtrials, 
    args.batchIndex, 
    args.tw_in_second))
print("\n...Done\n")

# if args.concat and args.batchIndex * args.batchNtrials == args.totalNtrials:
#     print("\n===== Concatenating {} tw{} =====\n".format(args.grb_name, args.tw_in_second))
#     import time
#     files = []
#     cnt = 0
#     do_concat = False
#     while True:
#         files = glob(output_folder + "{}_batchSize{}_batchIndex*_tw{}.pkl".format(
#             args.grb_name, 
#             args.batchNtrials, 
#             args.tw_in_second
#         ))
#         if len(files) == args.batchIndex + 1:
#             do_concat = True
#             break
#         time.sleep(180)
#         cnt += 1
#         if cnt > 100: 
#             break
#     if do_concat:
#         files = sorted(files)
#         df_list = [pd.read_pickle(x) for x in files]
#         out_df = pd.concat(df_list)
#         sout_df = out_df.astype(pd.SparseDtype("float", 0.0))
#         try:
#             sout_df.to_pickle(output_folder+"/{}_totalSize{}_tw{}.pkl".format(args.grb_name, 
#                                                                                 args.totalNtrials, 
#                                                                                 args.tw_in_second))
#             os.system("rm {}/*_batchSize*.pkl".format(output_folder))
#         except:
#             raise Exception("cannot save sparse concatenated pickle or cannot rm batch pickles.")
#     print("\n...Done\n")

