#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/env-shell.sh
#!/usr/bin/env python
r"""
Perform all-sky scans with no spatial prior information.
Record the TS (and number of true and fitted events around the
best fit location) for each pixel, creating a background-only TSD
at each point in the sky.
Combine this map with GBM healpix maps or grab a single pixel
for non-GBM bursts.
"""

import sys
# Put csky in my python path when using OSG
modules_dir = '/cvmfs/icecube.opensciencegrid.org/users/cjchen'
sys.path.append(modules_dir+'/csky/')
sys.path.append(modules_dir+'/greco_grb/')


import os
import sys
import numpy as np
import healpy as hp
import histlite as hl
import csky as cy
import pandas as pd
from scipy import sparse

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
%matplotlib inline
# %matplotlib notebook
from glob import glob
timer = cy.timing.Timer()
time = timer.time

###### Local Import ######
import SETTING
paths = SETTING.PATH()
print(paths)
USER = paths.USER
ICDATA_DIR = paths.ICDATA_DIR
DATA_DIR = paths.DATA_DIR
ANA_DIR = paths.ANA_DIR

from utils import *

print("Print statement after the imports")

