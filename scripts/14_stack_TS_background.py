"""
For the comparison plot: to stack TS, we need sob and n_b
"""

import pandas as pd
import sys
sys.path.append('../../')
import csky as cy
from greco_grb.scripts import SETTING
paths = SETTING.PATH()
print(paths)
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR
from greco_grb.scripts.utils import *
from glob import glob
import histlite as hl
import scipy

import argparse
######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Signal Trials",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--grb_idx_start", default=0, type=int, help="GRB idx start: [start, end)")
p.add_argument("--grb_idx_end", default=1, type=int, help="GRB idx end: [start, end)")
p.add_argument("--batchNtrials", default=1000, type=int, help="Number of trials")
p.add_argument("--batchIndex", default=0, type=int, help="Current batchIdx for this n_inj with this tw_in_second")
p.add_argument("--ncpu", default=4, type=int, help="Number of cores used")
args = p.parse_args()
###########################################################################
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
    'cut_n_sigma': 3,
}

cy.CONF.update(conf)
df = pd.read_pickle(DATA_DIR+"/grbwebgbm/grbweb_gbm_noHealpix_2268.pkl")




for grb_idx in range(args.grb_idx_start, args.grb_idx_end):
    grb_name = df.grb_name[grb_idx]
    grb_row = df.loc[df['grb_name'] == grb_name]
    ra = grb_row.ra
    dec = grb_row.dec
    t100 = grb_row.t100.values[0] * 86400.
    ### load grb healpix map
    healpix = np.load(DATA_DIR+"/grbwebgbm/healpix/{}_healpix_nside64.npy".format(grb_name))
    healpix = np.maximum(healpix,0)
    ########## healpix reduce (< instead of <=) ##########
    healpix[healpix < isf_healpix(healpix, q=0.99)] = 0
    healpix = healpix / np.sum(healpix)
    
    rng=np.random.default_rng(abs(java_hash(grb_name)))
    seeds = rng.integers(int(1e9), size=int(2e8))[args.batchNtrials*args.batchIndex: args.batchNtrials*(args.batchIndex + 1)]
    print("\n...Done seeding\n")
    ### for each tw_in_second
    for tw_idx, tw_in_second in enumerate([t100, 10,25,50,100,250,500]):
        timer = cy.timing.Timer()
        time = timer.time
        with time(f"Working on {grb_name} with TW{tw_in_second}..."):
            sob_list = []
            n_b_list =[]
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
            for no_trial in range(args.batchNtrials):
                sptr = cy.get_spatial_prior_trial_runner(conf=cy.CONF
                                             ,src_tr=src
                                             ,llh_priors=[healpix])

                ############## ref: osg_stacked_bg_tsd.py ########################
                ############## Get the pixelmax ##################################
                try:
                    _, ts, ns, ra_fit, dec_fit = sptr.get_one_fit(TRUTH=False, 
                                                                   mp_cpus=args.ncpu, 
                                                                   logging=False,
                                                                  seed=seeds[no_trial])
                except:
                    ts, ns, ra_fit, dec_fit = 0.0, 0.0, 0.0, 0.0
                if ns == 0.0:
                    ra_fit, dec_fit = grb_row.ra, grb_row.dec
                ############## Get SoB, n_b for point source: pixelmax ##########
                src_fit = cy.sources(
                                ra=ra_fit,
                                dec=dec_fit,
                                deg=True,
                                mjd=tw_start, 
                                sigma_t=np.zeros_like(tw), 
                                t_100=tw,  # in days
                                # prior=[hl.heal.HealHist(healpix)], # not important here
                                name=grb_name
                            )
                tr = cy.get_trial_runner(conf=cy.CONF, ana=ana, src=src_fit)
                L = tr.get_one_llh(TRUTH=False, seed=seeds[no_trial])
                try:
                    #### i=0 for bg, i=1 for sig; unblind so i=0.
                    SB_space = cy.inspect.get_space_eval(L, -1, 0)()[0]
                    ## unblinded, no need for (gamma=xxx)
                    SB_energy = cy.inspect.get_energy_eval(L, -1, 0)()[0]
                    SB = SB_space * SB_energy
                except:
                    SB = np.array([])
                n_b = L.llh_model.N_bg
                sob_list.append(SB)
                n_b_list.append(n_b)
        # Save sob_list and n_b_list for this GRB with this TW
        print(f"Saving SoB and n_b for {grb_name} with TW={tw_in_second:.2f}")
        output_path = cy.utils.ensure_dir(DATA_DIR+f"/csky_output/comparison/background/{grb_name}")
        np.save(f"{output_path}/{grb_name}_tw{tw_in_second:.2f}_batchIndex{args.batchIndex}_batchNtrials{args.batchNtrials}.npy", np.array(sob_list))
        np.save(f"{output_path}/{grb_name}_tw{tw_in_second:.2f}_batchIndex{args.batchIndex}_batchNtrials{args.batchNtrials}.npy", np.array(n_b_list))