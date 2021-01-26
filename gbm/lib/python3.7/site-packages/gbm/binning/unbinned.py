# unbinned.py: Module containing data binning functions for unbinned data
#
#     Authors: William Cleveland (USRA),
#              Adam Goldstein (USRA) and
#              Daniel Kocevski (NASA)
#
#     Portions of the code are Copyright 2020 William Cleveland and
#     Adam Goldstein, Universities Space Research Association
#     All rights reserved.
#
#     Written for the Fermi Gamma-ray Burst Monitor (Fermi-GBM)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import numpy as np


def bin_by_time(times, dt, tstart=None, tstop=None, time_ref=None):
    """Bins unbinned data to a specified temporal bin width.
    
    Args:
        times (np.array): The time of each event
        dt (float): The requested temporal bin width in seconds
        tstart (float, optional): The first bin edge time. Will use the first 
                                  event time if omitted.
        tstop: (float, optional): The last bin edge time. Will use the last 
                                  event time if omitted.
        time_ref (float, optional): 
            The reference time at which the binning will be based. If the set, 
            the binning will proceed starting at ``time_ref`` and moving forward 
            in time as well as starting at ``time_ref`` and moving backward in 
            time.  If not set, the binning will start at the beginning of the data.
    
    Returns:
        np.array: The edges of the binned data
    """
    assert dt > 0.0, "Requested bin width must be > 0.0 s"
    if time_ref is not None:
        time_ref = float(time_ref)

    if tstart is None:
        tstart = np.min(times)
    if tstop is None:
        tstop = np.max(times)

    # if we are using a reference time
    if time_ref is not None:
        pre_edges = np.arange(time_ref, tstart - dt, -dt)[::-1]
        post_edges = np.arange(time_ref, tstop + dt, dt)
        edges = np.concatenate((pre_edges[:-1], post_edges))
    else:
        edges = np.arange(tstart, tstop + dt, dt)

    return edges


def combine_into_one(times, tstart, tstop):
    """Bins unbinned data to a single bin.
        
    Args:
        times (np.array): The time of each event
        tstart (float): The first bin edge time. Will use the first 
                                  event time if omitted.
        tstop: (float): The last bin edge time. Will use the last 
                                  event time if omitted.
    
    Returns:
        np.array: The edges of the binned data
    """
    assert tstart < tstop, "The time range must be set in " \
                           "ascending order"

    # if no events, then the we have no counts
    if times.size == 0:
        return np.array([0]), np.array((tstart, tstop))

    if (np.min(times) > tstop) | (np.max(times) < tstart):
        raise ValueError("Requested time range is outside data range")

    time_range = (tstart, tstop)
    return np.array(time_range)


def combine_by_factor(times, old_edges, bin_factor, tstart=None, tstop=None):
    """Bins individual events to a multiple factor of bins given a 
    set of bin edges
    
    Args:
        times (np.array): The time of each event
        old_edges (np.array): The edges to be combined
        bin_factor (int): The number of bins to be combined
        tstart (float, optional): The first bin edge time. Will use the first 
                                  event time if omitted.
        tstop: (float, optional): The last bin edge time. Will use the last 
                                  event time if omitted.
    
    Returns:
        np.array: The edges of the binned data
    """

    assert bin_factor >= 1, "bin_factor must be a positive integer"
    bin_factor = int(bin_factor)

    # mask the old edges for the desired data range
    if tstart is None:
        tstart = old_edges[0]
    if tstop is None:
        tstop = old_edges[-1]
    old_edges = old_edges[old_edges >= tstart]
    old_edges = old_edges[old_edges <= tstop]

    # create the new edges
    new_edges = old_edges[::bin_factor]
    return new_edges


def bin_by_snr(times, back_rates, snr):
    """Bins unbinned data by SNR
        
    Args:
        times (np.array): The time of each event
        back_rates (np.array): The background rate at the time of each event
        snr (float): The signal-to-noise ratio
    
    Returns:
        np.array: The edges of the binned data
    """
    # get the background rates and differential counts at each event time
    back_counts = np.zeros_like(back_rates)
    back_counts[:-1] = back_rates[:-1] * (times[1:] - times[:-1])
    back_counts[-1] = back_counts[-2]

    num_events = len(times)
    istart = 0
    edges = []
    while True:
        # cumulative sum of the counts, background counts, and snr
        countscum = np.arange(1, num_events + 1 - istart)
        backgroundscum = np.cumsum(back_counts[istart:])
        snrcum = (countscum - backgroundscum) / np.sqrt(backgroundscum)
        # determine where to make the cut
        below_thresh = np.sum(snrcum <= snr)
        if below_thresh == 0:
            iend = istart
        else:
            iend = istart + below_thresh - 1
        edges.append(iend)
        if iend >= num_events - 1:
            break
        istart = iend + 1

    # get the finalized edges
    if edges[-1] != num_events - 1:
        edges.append(num_events - 1)
    edges = times[np.array(edges)]
    return edges


def time_to_spill(times, threshold):
    """Time-to-Spill Binning for an event list
    Bins an event list by accumulating counts until the set threshold is reached, 
    and then creating a new bin.

    Args:
        times (np.array): The time of each event
        threshold (int): The count threshold for the histogram
    
    Returns:
        np.array: The edges of the binned data
    """
    threshold = int(threshold)
    assert threshold > 0, "Threshold must be positive"
    # this is easy: take only the nth time (set by threshold) as the bin edge
    edges = times[::threshold]
    return edges


def bin_by_edges(times, time_edges):
    """Bins unbinned data by pre-defined edges. A rather trivial function :)
        
    Args:
        times (np.array): The time of each event
        time_edges (np.array): The pre-defined time edges
    
    Returns:
        np.array: The edges of the binned data
    """
    return time_edges
