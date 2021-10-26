"""
Some of my commonly used functions
"""
import numpy as np
from scipy import sparse

def java_hash(s):
    """
    A hash function implemented in Java style
    
    Parameters
    ----------
    s : str
        input string
    
    Returns
    -------
    x : scalar
        an integer (could be negative)
    
    See also
    -------
    http://hg.openjdk.java.net/jdk8/jdk8/jdk/file/687fd7c7986d/src/share/classes/java/lang/String.java#l1452
    """
    h = 0
    for c in s:
        h = (31 * h + ord(c)) & 0xFFFFFFFF
    return ((h + 0x80000000) & 0xFFFFFFFF) - 0x80000000


def isf_healpix(arr, q=0.999):
    """
    A customized isf for a healpix map (Inverse Survival Function (Inverse of SF))
    It is much safer than skylab.priors.pixel_cut().
    
    Parameters
    ----------
    arr : array_like
        A healpix map. normally we have np.sum(arr) = 1. but it's not neccessary.
        
    q : array_like (0.0 ~ 1.0, default=0.999)
        upper tail probability
    
    Returns
    -------
    x : scalar
        The value in arr corresponding to the upper tail probability q.
    
    See also
    -------
    https://docs.scipy.org/doc/scipy/reference/tutorial/stats/discrete.html#inverse-survival-function
    https://docs.scipy.org/doc/scipy/reference/tutorial/stats.html#common-methods
    """
    try:
        import numpy as np
    except:
        raise ImportError("Numpy cannot imported")
    if q == 1:
        return 0
    sorted_arr = np.sort(arr)
    cumsum = np.cumsum(sorted_arr)
    cdf = cumsum / cumsum[-1]
    # sf = 1 - cdf # no need
    # np.searchsorted(a,v): a has to be ascending
    idx = np.searchsorted(cdf, 1-q, side="left")
    # we return the nearest (to q) one
    if idx > 0 and (idx == len(sorted_arr) or abs(1-q - cdf[idx-1]) < abs(1-q - cdf[idx])):
        return sorted_arr[idx-1]
    return sorted_arr[idx]
    
def stay_awake(interval=180):
    """
    A stupid function that can print the running time after each interval, so that it keeps the notebook awake.
    
    Parameters
    ----------
    interval : scalar
        The unit is second.
    
    Returns
    -------
    N/A
    
    See also
    -------
    N/A
    
    """
    try:
        import time
    except:
        raise Exception("caanot import time")
    cnt = 0
    while True:
        print("{:d} min".format(int(cnt*interval/60.)), end = ' -> ')
        time.sleep(interval)
        cnt += 1
        
def ensure_dir(dirname):
    """
    Make sure ``dirname`` exists and is a directory.
    
    Parameters
    ----------
    dirname : str
        The pathname of the directory
    
    Returns
    -------
    s : str
        The pathname of the directory
    
    See also
    -------
    An implementation from csky: csky.utils.ensure_dir()
    """
    try:
        import os
        from os import errno
    except:
        import errno
    if not os.path.isdir(dirname):
        try:
            os.makedirs(dirname)   # throws if exists as file
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    return dirname


def getCombinedMaxTSfromFile(healpix, file):
    """
    Parameters
    ----------
    healpix : array_like
        A prior healpix map. shape [hp.nside2pix(64), ]
        
    file : str
        Path to the no_prior background trial file. One file could contain
        multiple trials
    
    Returns
    -------
    TSs : list
        A list of maximum TSs from the combined backgrounds.
    
    See also
    -------
    N/A
    """
    scans = sparse.load_npz(file)
    healpix = np.maximum(healpix, 1e-15)
    tss = scans.copy()
    # non-overlapping region will have negative TS anyways
    tss.data += 2. * (np.log(healpix[scans.indices]) - np.log(np.max(healpix)))
    tss.data[~np.isfinite(tss.data)] = 0
    tss.data[tss.data < 0] = 0
    TSs = tss.max(axis=1).A
    TSs = np.hstack(np.array(TSs))
    return TSs

