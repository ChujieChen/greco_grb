import os
import sys
import numpy as np
import histlite as hl
import csky as cy
import pandas as pd

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
p = argparse.ArgumentParser(description="Unbindling script",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_idx_start", default=0, type=int, help="GRB idx start: [start, end)")
p.add_argument("--grb_idx_end", default=1, type=int, help="GRB idx end: [start, end)")
p.add_argument("--ncpu", default=1, type=int, help="Number of CPUs")
p.add_argument("--out_dir", default="/storage/home/hhive1/cchen641/data/icecube/data/greco_grb/data/csky_output", type=str, help="Output root directory")
args = p.parse_args()

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
                      , load_sig=True)  # to save memory: use False

############# used for basic trial_runner
# conf = {
#     'ana': ana,
#     #### llh basics: csky.conf
#     'space': 'prior', # ps/fitps/template/prior
#     'time': 'transient', # utf/lc/transient
#     'energy': 'customflux', # fit/customflux
#     'flux': cy.hyp.PowerLawFlux(2.5),
#     #### inj.py - prior has some duplications against space's prior
#     'sig': 'transient', # ps/tw/lc/transient/template/prior
#     'full_sky': True,
#     'extended': True,
#     'cut_n_sigma': 3,
#     "TRUTH": True,
#     'mp_cpus': args.ncpu
# }

############# used for spatial_prior_trial_runner
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

df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")

# for each grb
for grb_idx in range(args.grb_idx_start, args.grb_idx_end):
    # load grb info
    grb_name = df.grb_name[grb_idx]
    print(f"Woring on {grb_name}")
    grb_row = df.loc[df['grb_name'] == grb_name]
    ra = grb_row.ra
    dec = grb_row.dec
    # load grb healpix map
    healpix = np.load(DATA_DIR+"/grbwebgbm/healpix/{}_healpix_nside64.npy".format(grb_name))
    healpix = np.maximum(healpix,0)
    ########## healpix reduce (< instead of <=) ##########
    healpix[healpix < isf_healpix(healpix, q=0.99)] = 0
    healpix = healpix / np.sum(healpix)
    # for each tw_in_second
    tss, nss = [], []
    for tw_in_second in [10,25,50,100,250,500]:
        tw = tw_in_second/86400.
        tw_start = grb_row.t_center - 0.5*tw
        src = cy.sources(
            ra=ra,
            dec=dec,
            deg=True,
            mjd=tw_start, 
            sigma_t=np.zeros_like(tw), 
            t_100=tw,  # in days
            prior=[hl.heal.HealHist(healpix)],
            name=grb_name
        )
        print(f"tw_in_second = {tw_in_second}")
        # tr = cy.get_trial_runner(conf=cy.CONF, ana=ana, src=src)
        # result = ts, ns = tr.get_one_fit(TRUTH=True)
        sptr = cy.get_spatial_prior_trial_runner(conf=cy.CONF
                                         ,src_tr=src
                                         ,llh_priors=[healpix])
        try:
            result = _, ts, ns, ra, dec = sptr.get_one_fit(TRUTH=True, mp_cpus=args.ncpu, logging=False)
        except:
            ts, ns = 0.0, 0.0
        tss.append(ts)
        nss.append(ns)
    # get pvals
    bg_files = sorted(glob(ANA_DIR + f"/allsky_scan/with_prior_background/tw*/{grb_name}_*.npz"), 
                 key=lambda x: int(x[x.find("/tw")+3:x.find(f"/{grb_name}")]))
    assert len(bg_files) > 0, print(f"cannot find bkg files of {grb_name}")
    bgs = np.array([sparse.load_npz(bg_file).toarray()[0] for bg_file in bg_files])
    bgs_sorted = np.apply_along_axis(sorted, 1, bgs)
    pvals = []
    for ts, bg_sorted in zip(tss, bgs_sorted):
        pvals.append((bg_sorted.size - np.searchsorted(bg_sorted, ts, side='left')) / bg_sorted.size)
    tss = np.array(tss)
    nss = np.array(nss)
    pvals = np.array(pvals)
    try:
        ts_ns_pval = np.array(
            list(map(tuple, np.transpose([tss,nss,pvals]))),
            dtype=np.dtype(
                [
                    ('ts', np.float32), 
                    ('ns', np.float32), 
                    ('pval', np.float32)
                ]
            )
        )
        output_folder = cy.utils.ensure_dir(args.out_dir+f"/unblind/ts_ns_p")
        print(f"Saving {output_folder}/{grb_name}_ts_ns_p.npy")
        np.save(output_folder + f"/{grb_name}_ts_ns_p.npy", ts_ns_pval)
    except:
        print("An exception occurs when saving.")
        
print("All done")