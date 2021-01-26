# unbinned.py: Module containing background fitting classes for unbinned data
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
from scipy.interpolate import interp1d


class NaivePoisson:
    """ A class to estimate the background of un-binned data using the naive Poisson 
    maximum likelihood.
    
    This method is approximately equivalent to sliding a window of fixed length 
    through un-binned data and calculating the Poisson maximum likelihood for 
    the rate. The rate estimate is applied to the center of the sliding window, 
    therefore, the amount of data equivalent to half of the sliding window at 
    the beginning and half of the window at the end of the data is a constant.
    
    Note:
        This naive approach assumes there is either no 'strong' signal in the \
         data, or the presence of a 'weaker' signal has a duration much less \
         than the window width of the sliding window.    
    
    Parameters:
        counts (list of np.array): A list of length numchans, and each element 
                                   of the list is an array of event times in 
                                   that channel.
    """

    def __init__(self, times):
        self._times = times
        self._numchans = len(times)
        self._min_dt = 2e-6
        self._window_width = None
        self._actual_widths = None
        self._rates = None

    def fit(self, window_width=100.0, fast=True):
        """Fit the data via Naive Poisson Maximum Likelihood
        
        Args:
            window_width (float): 
                The width of the sliding window in seconds. 
                
                Note: 
                    If the range of the data is shorter than ``window_width``,
                    the ``window_width`` will automatically be shortened to the 
                    range of the data.
            fast (bool): If True, then will use the fast approximation of the 
                     algorithm that allows the ``window_width`` to change 
                     throughout the data (the number of counts in the window 
                     is constant). If False, uses the exact algorithm with a 
                     fixed window, but is much slower.
        """
        self._window_width = window_width
        actual_widths = []
        rates = []
        uncerts = []
        for i in range(self._numchans):
            if fast:
                r, u, w = self._fit_one_fast(i)
            else:
                r, u, w = self._fit_one_exact(i)
            rates.append(r)
            uncerts.append(u)
            actual_widths.append(w)

        self._actual_widths = actual_widths
        self._rates = rates

    def _fit_one_exact(self, channel):
        """Fit a single channel of event data.  This is the exact (not
        approximate) algorithm that uses a fixed window duration throughout the
        data, except where the window must be necessarily truncated at the ends
        of the data.  This function is much slower than the approximate version.
        
        Args:
            channel (int): The energy channel
        
        Returns:        
            (np.array, np.array, np.array): The background rates, uncertainty, \
                                            and actual window width
        """
        window_width = self._window_width
        events = self._times[channel]

        num_events = len(events)
        if num_events == 0:
            return (np.array([]), np.array([]), np.array([]))

        # get the different parts of the array
        # pre: before the full sliding window width
        # mid: during the full sliding window width
        # post: after the full sliding window width
        pre_mask = (events < events[0] + window_width / 2.0)
        post_mask = (events > events[-1] - window_width / 2.0)
        mid_mask = (~pre_mask & ~post_mask)

        # do the sliding window
        idx = np.sum(pre_mask)
        mid_bins = [
            np.sum(np.abs(events - events[i + idx]) <= window_width / 2.0) \
            for i in range(np.sum(mid_mask))]
        mid_rates = np.array(mid_bins) / window_width
        mid_uncert = np.sqrt(mid_rates / window_width)

        # now considering the pre- and post-full-window-width data:
        # assume a constant rate at either end, but shorten the window
        # appropriately

        pre_events = events[pre_mask]
        pre_dt = (pre_events - events[0])
        pre_rates = np.full(np.sum(pre_mask) - 1, mid_rates[0])
        pre_uncert = np.sqrt(pre_rates / pre_dt[1:])

        idx = num_events - np.sum(post_mask)
        post_events = events[post_mask]
        post_dt = events[-1] - post_events
        post_rates = np.full(np.sum(post_mask) - 1, mid_rates[-1])
        post_uncert = np.sqrt(post_rates / post_dt[:-1])

        # put it all together
        brates = np.hstack((pre_rates, mid_rates, post_rates))
        uncert = np.hstack((pre_uncert, mid_uncert, post_uncert))
        dt = np.hstack((pre_dt[1:], np.full(np.sum(mid_mask), window_width),
                        post_dt[:-1]))

        return brates, uncert, dt

    def _fit_one_fast(self, channel):
        """Fit a single channel of event data.  This performs an approximation
        to the NaivePoisson algorithm to considerably speed up the computation.
        Instead of assuming a fixed window duration throughout the data, the 
        window is initialized by determining the number of counts in the first
        window, and fixing the window width by the total number of counts in
        the data.  This allows the window duration to change, and is similar
        to smoothing the data.  For slowly varying data, this is a good 
        approximation. 
        
        Args:
            channel (int): The energy channel
        
        Returns:        
            (np.array, np.array, np.array): The background rates, uncertainty, \
                                            and actual window width
        """
        window_width = self._window_width
        events = self._times[channel]
        num_events = len(events)
        if num_events == 0:
            return (np.array([]), np.array([]), np.array([]))

        # get extent of window at the start of the data
        end_idx = np.sum(events <= events[0] + window_width)
        # from the middle of the sliding window
        mid_idx = int(np.floor(end_idx / 2.0))

        # number of 'bins' containing a count (= number of counts in window)
        full_bins = end_idx
        # number of 'bins' containing a count not in the window
        num_back_bins = num_events - end_idx
        array_back = np.arange(num_back_bins)

        # actual window widths and the rate
        dts = (events[end_idx + array_back] - events[array_back])
        back_rates = full_bins / dts

        # Now we want to estimate the rate at the ends. This means that the 
        # sliding window will truncate at either end until it's the width of 
        # one event.  The main concern is to calculate the width of the 
        # truncated window correctly so that the uncertainty in the rate will 
        # properly increase as the window width decreases

        pre_full_bins = np.arange(mid_idx) + 1
        pre_dts = events[2 * (pre_full_bins[1:] - 1)] - events[0]
        # pre_back_rates = np.full(mid_idx-1, back_rates[0])
        pre_back_rates = (2.0 * pre_full_bins[1:]) / pre_dts

        post_full_bins = pre_full_bins[::-1]
        post_dts = events[-1] - events[-1 - 2 * (post_full_bins[:-1] - 1)]
        post_back_rates = (2.0 * post_full_bins[:-1]) / post_dts
        # post_back_rates = np.full(mid_idx-1, back_rates[-1])

        # put all of it together
        brates = np.hstack((pre_back_rates, back_rates, post_back_rates))
        dts = np.hstack((pre_dts, dts, post_dts))

        # this is if we ended up with an odd number of events
        if num_events - 2 > brates.shape[0]:
            brates = np.append(brates, brates[-1])
            dts = np.append(dts, dts[-1])

        # now the uncertainty
        uncert = np.sqrt(brates / dts)

        return (brates, uncert, dts)

    def interpolate(self, tstart, tstop):
        """Interpolate the background at the given times
        
        Args:
            tstart (np.array): The start times of the bins to interpolate
            tstop (np.array): The end times of the bins to interpolate
        
        Returns:        
            (np.array, np.array): The interpolated model value and model \
                                  uncertainty in each bin
        """
        times = (tstart + tstop) / 2.0
        rates = []
        uncert = []
        for i in range(self._numchans):
            rates_interp = interp1d(self._times[i][1:-1], self._rates[i],
                                    fill_value='extrapolate')
            width_interp = interp1d(self._times[i][1:-1],
                                    self._actual_widths[i],
                                    fill_value='extrapolate')
            r = rates_interp(times)
            widths = width_interp(times)
            uncert.append(np.sqrt(r / widths))
            rates.append(r)

        rates = np.array(rates).T
        uncert = np.array(uncert).T
        return rates, uncert
