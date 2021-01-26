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
import unittest
import numpy as np
from gbm.data.primitives import EventList, TimeBins, EnergyBins, TimeEnergyBins
from gbm.data.primitives import TimeRange, GTI
from gbm.binning.binned import combine_by_factor
from gbm.binning.unbinned import bin_by_time

class TestEventList(unittest.TestCase):
    events = [(1.4, 1), (1.5, 0), (2.3, 0), (2.4, 1), (3.9, 0),
              (6.8, 0), (7.6, 0), (7.7, 2), (8.6, 1), (9.8, 3)]
    events = np.array(events, dtype=[('TIME', '>f8'), ('PHA', '>i2')])
    ebounds = [(0, 4.2, 5.2), (1, 5.2, 6.1), (2, 6.1, 7.0), (3, 7.0, 7.8)]
    ebounds = np.array(ebounds, dtype=[('CHANNEL', '>i2'), ('E_MIN', '>f4'), ('E_MAX', '>f4')])
    
    def test_attributes(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        self.assertEqual(el.size, 10)
        self.assertCountEqual(el.time, self.events['TIME'])
        self.assertCountEqual(el.pha, self.events['PHA'])
        self.assertCountEqual(el.emin, self.ebounds['E_MIN'])
        self.assertCountEqual(el.emax, self.ebounds['E_MAX'])
        self.assertEqual(el.numchans, 4)
        self.assertEqual(el.time_range, (self.events['TIME'][0], self.events['TIME'][-1])) 
        self.assertEqual(el.channel_range, (self.ebounds['CHANNEL'][0], self.ebounds['CHANNEL'][-1])) 
        self.assertEqual(el.energy_range, (self.ebounds['E_MIN'][0], self.ebounds['E_MAX'][-1])) 
        self.assertAlmostEqual(el.get_exposure(), 8.39997, places=4)
        self.assertIsInstance(el.count_spectrum, EnergyBins)
    
    def test_channel_slice(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el2 = el.channel_slice(0, 1)
        self.assertEqual(el2.size, 8)
        self.assertEqual(el2.time_range, (1.4, 8.6))

    def test_count_spectrum(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        bins = el.count_spectrum
        self.assertEqual(bins.size, 4)
        self.assertCountEqual(bins.counts, np.array([5,3,1,1]))
    
    def test_energy_slice(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el2 = el.energy_slice(5.0, 6.2)
        self.assertEqual(el2.size, 9)
        self.assertEqual(el2.time_range, (1.4, 8.6))
        
    def test_time_slice(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el2 = el.time_slice(1.0, 3.0)
        self.assertEqual(el2.size, 4)
        self.assertEqual(el2.time_range, (1.4, 2.4))
    
    def test_sort(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el.sort('PHA')
        self.assertEqual(el.pha[0], 0)
    
    def test_bin_time(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        bins1 = el.bin(bin_by_time, 1.0)
        bins2 = el.bin(bin_by_time, 1.0, tstart=2.0, tstop=8.0)
        self.assertEqual(bins1.size, (9, 4))
        self.assertEqual(bins2.size, (6, 4))

    def test_rebin_energy(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el2 = el.rebin_energy(combine_by_factor, 2)
        self.assertEqual(el2.numchans, 2)
    
    def test_merge(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        el1 = el.time_slice(1.0, 5.0)
        el2 = el.time_slice(5.0, 10.0)
        el3 = EventList.merge([el1, el2], sort_attrib='TIME')
        self.assertEqual(el3.size, el.size)
        self.assertCountEqual(el3.time, el.time)
    
    def test_exposure(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        exposure = el.get_exposure()
        self.assertAlmostEqual(exposure, (9.8-1.4)-(10.0*2.6e-6), places=5)
        exposure = el.get_exposure(time_ranges=(0.0, 5.0))
        self.assertAlmostEqual(exposure, (5.0)-(5.0*2.6e-6), places=5)
        exposure = el.get_exposure(time_ranges=[(0.0, 2.0), (5.0, 9.0)])
        self.assertAlmostEqual(exposure, (2.0+4.0)-(6.0*2.6e-6))
    
    def test_from_lists(self):
        times = [1.4, 1.5, 2.3, 2.4, 3.9, 6.8, 7.6, 7.7, 8.6, 9.8]
        phas = [1,0,0,1,0,0,0,2,1,3]
        chan_lo = [4.2, 5.2, 6.1, 7.0]
        chan_hi = [5.2, 6.1, 7.0, 7.8]
        el = EventList.from_lists(times, phas, chan_lo, chan_hi)
        el2 = EventList.from_fits_array(self.events, self.ebounds)
        self.assertCountEqual(el.time, el2.time)
        self.assertCountEqual(el.pha, el2.pha)
        
    def test_errors(self):
        el = EventList.from_fits_array(self.events, self.ebounds)
        with self.assertRaises(ValueError):
            el.sort('PIGLET')
        

class TestTimeBins(unittest.TestCase):
    counts = np.array([113, 94, 103, 100, 98, 115, 101, 86, 104, 102])
    tstart = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    tstop = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    exposure = np.array([0.99971, 0.99976, 0.99973, 0.99974, 0.99975, 
                         0.99970, 0.99973, 0.99978, 0.99973, 0.99973])

    def test_attributes(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        self.assertCountEqual(bins.centroids, (self.tstart+self.tstop)/2.0)
        self.assertCountEqual(bins.counts, self.counts)
        self.assertCountEqual(bins.count_uncertainty, np.sqrt(self.counts))
        self.assertCountEqual(bins.exposure, self.exposure)
        self.assertCountEqual(bins.hi_edges, self.tstop)
        self.assertCountEqual(bins.lo_edges, self.tstart)
        self.assertEqual(bins.range, (0.0, 10.0))
        self.assertCountEqual(bins.rates, self.counts/self.exposure)
        self.assertCountEqual(bins.rate_uncertainty, np.sqrt(self.counts)/self.exposure)
        self.assertEqual(bins.size, 10)
        self.assertCountEqual(bins.widths, self.tstop-self.tstart)
    
    def test_closest_edge(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        self.assertEqual(bins.closest_edge(1.3), self.tstart[1])        
        self.assertEqual(bins.closest_edge(1.3, which='high'), self.tstart[2])        
        self.assertEqual(bins.closest_edge(1.7, which='low'), self.tstart[1])
        
        self.assertEqual(bins.closest_edge(-1.0), self.tstart[0])        
        self.assertEqual(bins.closest_edge(-1.0, which='low'), self.tstart[0])        
        self.assertEqual(bins.closest_edge(-1.0, which='high'), self.tstart[0])        

        self.assertEqual(bins.closest_edge(11.0), self.tstop[-1])        
        self.assertEqual(bins.closest_edge(11.0, which='low'), self.tstop[-1])        
        self.assertEqual(bins.closest_edge(11.0, which='high'), self.tstop[-1])        
    
    def test_slice(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins2 = bins.slice(2.5, 7.5)
        self.assertEqual(bins2.size, 6)
        self.assertCountEqual(bins2.counts, self.counts[2:8])
        self.assertEqual(bins2.range, (2.0, 8.0))
        
    def test_rebin(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins2 = bins.rebin(combine_by_factor, 2)
        bins3 = bins.rebin(combine_by_factor, 2, tstart=2.5, tstop=7.5)
        self.assertEqual(bins2.size, 5)
        self.assertCountEqual(bins2.counts, np.sum(self.counts.reshape(-1,2), axis=1))
        self.assertCountEqual(bins2.exposure, np.sum(self.exposure.reshape(-1,2), axis=1))
        self.assertEqual(bins3.size, 8)
    
    def test_merge(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins1 = bins.slice(0.1, 4.9)
        bins2 = bins.slice(5.1, 9.9)
        bins3 = TimeBins.merge([bins1, bins2])
        self.assertCountEqual(bins3.counts, bins.counts)
        self.assertCountEqual(bins3.exposure, bins.exposure)
        self.assertCountEqual(bins3.lo_edges, bins.lo_edges)
        self.assertCountEqual(bins3.hi_edges, bins.hi_edges)
    
    def test_sum(self):
        bins1 = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins2 = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins3 = TimeBins.sum([bins1, bins2])
        self.assertCountEqual(bins3.counts, self.counts*2)
        self.assertCountEqual(bins3.exposure, self.exposure)
    
    def test_contiguous(self):
        tstart = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.5, 6.0, 7.0, 8.0, 9.0])
        exposure = np.copy(self.exposure)
        exposure[5] = 0.4998
        bins = TimeBins(self.counts, tstart, self.tstop, exposure)
        split_bins = bins.contiguous_bins()
        self.assertEqual(len(split_bins), 2)
        self.assertEqual(split_bins[0].range, (0.0, 5.0))
        self.assertEqual(split_bins[1].range, (5.5, 10.0))

    def test_errors(self):
        bins = TimeBins(self.counts, self.tstart, self.tstop, self.exposure)
        bins2 = TimeBins(self.counts[1:], self.tstart[1:], self.tstop[1:], self.exposure[1:])
        with self.assertRaises(TypeError):
            bins.exposure = 1.0
        with self.assertRaises(TypeError):
            bins.counts = 1.0
        with self.assertRaises(TypeError):
            bins.lo_edges = 1.0
        with self.assertRaises(TypeError):
            bins.hi_edges = 1.0
        with self.assertRaises(AssertionError):
            TimeBins.sum([bins, bins2])


class TestEnergyBins(unittest.TestCase):
    counts = np.array([113, 94, 103, 100, 98, 115, 101, 86, 104, 102])
    emin = np.array([4.2, 5.2, 6.1, 7.0, 7.8, 8.7, 9.5, 10.4, 11.4, 12.4])
    emax = np.array([5.2, 6.1, 7.0, 7.8, 8.7, 9.5, 10.4, 11.4, 12.4, 13.5])
    exposure = np.full(10, 9.9)
    
    def test_attributes(self):
        bins = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        self.assertCountEqual(bins.centroids, np.sqrt(self.emin*self.emax))
        self.assertCountEqual(bins.counts, self.counts)
        self.assertCountEqual(bins.count_uncertainty, np.sqrt(self.counts))
        self.assertCountEqual(bins.exposure, self.exposure)
        self.assertCountEqual(bins.hi_edges, self.emax)
        self.assertCountEqual(bins.lo_edges, self.emin)
        self.assertEqual(bins.range, (4.2, 13.5))
        self.assertCountEqual(bins.rates, self.counts/(self.exposure*(self.emax-self.emin)))
        self.assertCountEqual(bins.rate_uncertainty, np.sqrt(self.counts)/(self.exposure*(self.emax-self.emin)))
        self.assertEqual(bins.size, 10)
        self.assertCountEqual(bins.widths, self.emax-self.emin)
    
    def test_slice(self):
        bins = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins2 = bins.slice(5.5, 10.5)
        self.assertEqual(bins2.size, 7)
        self.assertCountEqual(bins2.counts, self.counts[1:8])
        self.assertEqual(bins2.range, (5.2, 11.4))
    
    def test_rebin(self):
        bins = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins2 = bins.rebin(combine_by_factor, 2)
        bins3 = bins.rebin(combine_by_factor, 2, emin=5.5, emax=9.7)
        self.assertEqual(bins2.size, 5)
        self.assertCountEqual(bins2.counts, np.sum(self.counts.reshape(-1,2), axis=1))
        self.assertCountEqual(bins2.exposure, self.exposure[:5])
        self.assertEqual(bins3.size, 8)
    
    def test_merge(self):
        bins = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins1 = bins.slice(4.5, 8.0)
        bins2 = bins.slice(9.0, 13.4)
        bins3 = EnergyBins.merge([bins1, bins2])
        self.assertCountEqual(bins3.counts, bins.counts)
        self.assertCountEqual(bins3.exposure, bins.exposure)
        self.assertCountEqual(bins3.lo_edges, bins.lo_edges)
        self.assertCountEqual(bins3.hi_edges, bins.hi_edges)
    
    def test_sum(self):
        bins1 = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins2 = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins3 = EnergyBins.sum([bins1, bins2])
        self.assertCountEqual(bins3.counts, self.counts*2)
        self.assertCountEqual(bins3.exposure, self.exposure*2)
    
    def test_errors(self):
        bins = EnergyBins(self.counts, self.emin, self.emax, self.exposure)
        bins2 = EnergyBins(self.counts[1:], self.emin[1:], self.emax[1:], self.exposure[1:])
        with self.assertRaises(TypeError):
            bins.exposure = 1.0
        with self.assertRaises(TypeError):
            bins.counts = 1.0
        with self.assertRaises(TypeError):
            bins.lo_edges = 1.0
        with self.assertRaises(TypeError):
            bins.hi_edges = 1.0
        with self.assertRaises(AssertionError):
            EnergyBins.sum([bins, bins2])
    

class TestTimeEnergyBins(unittest.TestCase):
    counts = np.array([[113, 94, 103], [100, 98, 115], [101, 86, 104]])
    tstart = np.array([1.0, 2.0, 3.0])
    tstop = np.array([2.0, 3.0, 4.0])
    exposure = np.array([0.99971, 0.99976, 0.99973])
    emin = np.array([4.2, 5.2, 6.1])
    emax = np.array([5.2, 6.1, 7.0])
    chan_widths = (emax-emin)
    
    def test_attributes(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        self.assertCountEqual(bins.chan_widths, self.emax-self.emin)
        self.assertCountEqual(bins.emax, self.emax)
        self.assertCountEqual(bins.emin, self.emin)
        self.assertCountEqual(bins.energy_centroids, np.sqrt(self.emax*self.emin))
        self.assertEqual(bins.energy_range, (4.2, 7.0))
        self.assertCountEqual(bins.exposure, self.exposure)
        self.assertEqual(bins.numchans, 3)
        self.assertEqual(bins.numtimes, 3)
        self.assertEqual(bins.size, (3,3))
        self.assertCountEqual(bins.time_centroids, (self.tstart+self.tstop)/2.0)
        self.assertEqual(bins.time_range, (1.0, 4.0))
        self.assertCountEqual(bins.time_widths, self.tstop-self.tstart)
        self.assertCountEqual(bins.tstart, self.tstart)
        self.assertCountEqual(bins.tstop, self.tstop)
        for i in range(3):
            self.assertCountEqual(bins.counts[:,i], self.counts[:,i])
            self.assertCountEqual(bins.count_uncertainty[:,i], np.sqrt(self.counts[:,i]))
            self.assertCountEqual(bins.rates[:,i], 
                      self.counts[:,i]/self.exposure)
            self.assertCountEqual(bins.rate_uncertainty[:,i],
                      np.sqrt(self.counts[:,i])/self.exposure)

    def test_integrate_energy(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        lc1 = bins.integrate_energy()
        lc2 = bins.integrate_energy(emin=4.5, emax=5.5)
        self.assertCountEqual(lc1.counts, np.sum(self.counts, axis=1))
        self.assertCountEqual(lc2.counts, np.sum(self.counts[:,:-1], axis=1))
        
    def test_integrate_time(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        spec1 = bins.integrate_time()
        spec2 = bins.integrate_time(tstart=1.5, tstop=2.5)
        self.assertCountEqual(spec1.counts, np.sum(self.counts, axis=0))
        self.assertCountEqual(spec2.counts, np.sum(self.counts[:-1,:], axis=0))

    def test_rebin_energy(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins2 = bins.rebin_energy(combine_by_factor, 2)
        bins3 = bins.rebin_energy(combine_by_factor, 2, emin=4.5, emax=6.0)
        self.assertEqual(bins2.numchans, 1)
        self.assertEqual(bins3.numchans, 3)

    def test_rebin_time(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins2 = bins.rebin_time(combine_by_factor, 2)
        bins3 = bins.rebin_time(combine_by_factor, 2, tstart=1.5, tstop=2.5)
        self.assertEqual(bins2.numtimes, 1)
        self.assertEqual(bins3.numtimes, 3)

    def test_slice_energy(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins2 = bins.slice_energy(4.5, 5.5)
        self.assertEqual(bins2.numchans, 2)
        self.assertCountEqual(bins2.counts.flatten(), self.counts[:,:2].flatten())
        self.assertEqual(bins2.energy_range, (4.2, 6.1))

    def test_slice_time(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins2 = bins.slice_time(1.5, 2.5)
        self.assertEqual(bins2.numtimes, 2)
        self.assertCountEqual(bins2.counts.flatten(), self.counts[:2,:].flatten())
        self.assertEqual(bins2.time_range, (1.0, 3.0))

    def test_closest_time(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        self.assertEqual(bins.closest_time_edge(2.3), 2.0)
        self.assertEqual(bins.closest_time_edge(2.3, which='low'), 2.0)
        self.assertEqual(bins.closest_time_edge(2.3, which='high'), 3.0)

    def test_closest_energy(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        self.assertEqual(bins.closest_energy_edge(5.7), 6.1)
        self.assertEqual(bins.closest_energy_edge(5.7, which='low'), 5.2)
        self.assertEqual(bins.closest_energy_edge(5.7, which='high'), 6.1)

    def test_contiguous_time(self):
        tstart = np.array([1.0, 2.5, 3.0])
        exposure = np.copy(self.exposure)
        exposure[1] = 0.4998
        bins = TimeEnergyBins(self.counts, tstart, self.tstop, exposure,
                              self.emin, self.emax)
        split_bins = bins.contiguous_time_bins()
        self.assertEqual(len(split_bins), 2)
        self.assertEqual(split_bins[0].time_range, (1.0, 2.0))
        self.assertEqual(split_bins[1].time_range, (2.5, 4.0))

    def test_contiguous_energy(self):
        emin = np.array([4.2, 5.7, 6.1])
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              emin, self.emax)
        split_bins = bins.contiguous_energy_bins()
        self.assertEqual(len(split_bins), 2)
        self.assertEqual(split_bins[0].energy_range, (4.2, 5.2))
        self.assertEqual(split_bins[1].energy_range, (5.7, 7.0))
    
    def test_exposure(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        self.assertEqual(bins.get_exposure(), self.exposure.sum())        
        self.assertEqual(bins.get_exposure(time_ranges=(1.5, 2.5)), 
                         self.exposure[0:2].sum())        
    
    def test_merge_energy(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins1 = bins.slice_energy(4.2, 5.5)
        bins2 = bins.slice_energy(6.0, 7.0)
        bins3 = TimeEnergyBins.merge_energy([bins1, bins2])
        self.assertCountEqual(bins3.counts.flatten(), bins.counts.flatten())
        self.assertCountEqual(bins3.exposure, bins.exposure)
        self.assertCountEqual(bins3.tstart, bins.tstart)
        self.assertCountEqual(bins3.tstop, bins.tstop)
        self.assertCountEqual(bins3.emin, bins.emin)
        self.assertCountEqual(bins3.emax, bins.emax)

    def test_merge_time(self):
        bins = TimeEnergyBins(self.counts, self.tstart, self.tstop, self.exposure,
                              self.emin, self.emax)
        bins1 = bins.slice_time(1.5, 2.5)
        bins2 = bins.slice_time(3.1, 4.0)
        bins3 = TimeEnergyBins.merge_time([bins1, bins2])
        self.assertCountEqual(bins3.counts.flatten(), bins.counts.flatten())
        self.assertCountEqual(bins3.exposure, bins.exposure)
        self.assertCountEqual(bins3.tstart, bins.tstart)
        self.assertCountEqual(bins3.tstop, bins.tstop)
        self.assertCountEqual(bins3.emin, bins.emin)
        self.assertCountEqual(bins3.emax, bins.emax)


class TestTimeRange(unittest.TestCase):
    tstart = -5.0
    tstop = 10.0
    time_range = TimeRange(tstart, tstop)
    
    def test_attributes(self):
        self.assertEqual(self.time_range.tstart, self.tstart)
        self.assertEqual(self.time_range.tstop, self.tstop)
        self.assertEqual(self.time_range.duration, 15.0)
        self.assertEqual(self.time_range.center, 2.5)
        tstart, tstop = self.time_range.as_tuple()
        self.assertEqual(tstart, self.tstart)
        self.assertEqual(tstop, self.tstop)
    
    def test_contains(self):
        c = self.time_range.contains(1.0)
        self.assertTrue(c)
        c = self.time_range.contains(-10.0)
        self.assertFalse(c)
        c = self.time_range.contains(10.0)
        self.assertTrue(c)
        c = self.time_range.contains(10.0, inclusive=False)
        self.assertFalse(c)
        c = self.time_range.contains(1.0, inclusive=True)
        self.assertTrue(c)
        c = self.time_range.contains(-10.0, inclusive=True)
        self.assertFalse(c)
    
    def test_union(self):
        tr = TimeRange.union(self.time_range, TimeRange(5.0, 20.0))
        self.assertEqual(tr.tstart, -5.0)
        self.assertEqual(tr.tstop, 20.0)
        tr = TimeRange.union(self.time_range, TimeRange(-10.0, 0.0))
        self.assertEqual(tr.tstart, -10.0)
        self.assertEqual(tr.tstop, 10.0)
        tr = TimeRange.union(self.time_range, TimeRange(100, 200.0))
        self.assertEqual(tr.tstart, -5.0)
        self.assertEqual(tr.tstop, 200.0)
        tr = TimeRange.union(self.time_range, TimeRange(0.0, 5.0))
        self.assertEqual(tr.tstart, -5.0)
        self.assertEqual(tr.tstop, 10.0)
        tr = TimeRange.union(self.time_range, TimeRange(-100.0, 100.0))
        self.assertEqual(tr.tstart, -100.0)
        self.assertEqual(tr.tstop, 100.0)

    def test_intersection(self):
        tr = TimeRange.intersection(self.time_range, TimeRange(5.0, 20.0))
        self.assertEqual(tr.tstart, 5.0)        
        self.assertEqual(tr.tstop, 10.0)        
        tr = TimeRange.intersection(self.time_range, TimeRange(-10.0, 0.0))
        self.assertEqual(tr.tstart, -5.0)        
        self.assertEqual(tr.tstop, 0.0)        
        tr = TimeRange.intersection(self.time_range, TimeRange(1.0, 5.0))
        self.assertEqual(tr.tstart, 1.0)
        self.assertEqual(tr.tstop, 5.0)
        tr = TimeRange.intersection(self.time_range, TimeRange(-100.0, 100.0))
        self.assertEqual(tr.tstart, -5.0)
        self.assertEqual(tr.tstop, 10.0)
        tr = TimeRange.intersection(self.time_range, TimeRange(50.0, 100.0))
        self.assertIsNone(tr)

    def test_errors(self):
        with self.assertRaises(TypeError):
            self.time_range.contains('')
        with self.assertRaises(TypeError):
            tr = TimeRange(1.0, 'pennywise')
    

class TestGti(unittest.TestCase):
    time_list = [(1.0, 5.0), (20.0, 30.0), (50.0, 100.0)]
    gti = GTI.from_list(time_list)

    def test_attributes(self):
        self.assertEqual(self.gti.num_intervals, 3)
        tstart, tstop = self.gti.range
        self.assertEqual(tstart, 1.0)
        self.assertEqual(tstop, 100.0)
    
    def test_as_list(self, gti=None, time_list=None):
        if gti is None:
            gti = self.gti
        if time_list is None:
            time_list = self.time_list
        self.assertCountEqual(gti.as_list(), time_list)
    
    def test_insert(self):
        gti2 = GTI.from_list(self.time_list[:-1])
        gti2.insert(*self.time_list[-1])
        self.test_as_list(gti=gti2)
        gti2 = GTI.from_list(self.time_list[1:])
        gti2.insert(*self.time_list[0])
        self.test_as_list(gti=gti2)
        gti2 = GTI.from_list([(1.0, 5.0), (50.0, 100.0)])
        gti2.insert(*self.time_list[1])
        self.test_as_list(gti=gti2)
        
        gti2 = GTI.from_list(self.time_list)
        gti2.insert(2.0, 15.0)
        self.test_as_list(gti=gti2, 
                          time_list=[(1.0, 15.0), (20.0, 30.0), (50.0, 100.0)])
        gti2 = GTI.from_list(self.time_list)
        gti2.insert(2.0, 25.0)
        self.test_as_list(gti=gti2, time_list=[(1.0, 30.0), (50.0, 100.0)])
        gti2 = GTI.from_list(self.time_list)
        gti2.insert(-10.0, 2.0)
        self.test_as_list(gti=gti2, 
                          time_list=[(-10.0, 5.0), (20.0, 30.0), (50.0, 100.0)])
    
    def test_merge(self):
        gti2 = GTI.from_list([(8.0, 15.0), (150.0, 200.0)])
        gti3 = GTI.merge(self.gti, gti2)
        self.test_as_list(gti=gti3, time_list=[(1.0, 5.0), (8.0, 15.0), 
                                               (20.0, 30.0), (50.0, 100.0),
                                               (150.0, 200.0)])
        
        gti2 = GTI.from_list([(2.0, 15.0), (25.0, 75.0)])
        gti3 = GTI.merge(self.gti, gti2)
        self.test_as_list(gti=gti3, time_list=[(1.0, 15.0), (20.0, 100.0)])

        gti2 = GTI.from_list([(-10.0, 2.0), (75.0, 90.0)])
        gti3 = GTI.merge(self.gti, gti2)
        self.test_as_list(gti=gti3, time_list=[(-10.0, 5.0), (20.0, 30.0),
                                               (50.0, 100.0)])

    def test_contains(self):
        self.assertTrue(self.gti.contains(2.0))
        self.assertFalse(self.gti.contains(10.0))
        self.assertTrue(self.gti.contains(5.0))
        self.assertFalse(self.gti.contains(5.0, inclusive=False))

    def test_from_boolean_mask(self):
        times = np.arange(10)+1
        mask = np.array([0,0,0,1,1,1,0,0,0,0], dtype=bool)
        gti = GTI.from_boolean_mask(times, mask)
        self.assertEqual(gti.as_list()[0], (1.0, 3.0))
        self.assertEqual(gti.as_list()[1], (7.0, 10.0))
    
    def test_errors(self):
        with self.assertRaises(TypeError):
            self.gti.contains('')
        with self.assertRaises(TypeError):
            self.gti.insert(1.0, 'fivethirtyeight')

if __name__ == '__main__':
    unittest.main()
