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
p.add_argument("--out_dir", default="/storage/home/hhive1/cchen641/data/icecube/data/greco_grb/data/csky_output", type=str, help="Output root directory")
args = p.parse_args()
