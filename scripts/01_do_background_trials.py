import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import healpy as hp
import histlite as hl
import csky as cy
import pandas as pd

# %matplotlib inline
# %matplotlib notebook

from glob import glob

timer = cy.timing.Timer()
time = timer.time