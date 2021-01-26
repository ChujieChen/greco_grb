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
import os
from unittest import TestCase

from gbm.data.phaii import TTE
from gbm.binning.binned import combine_by_factor
from gbm.binning.unbinned import bin_by_time

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestTTE(TestCase):
    filename = os.path.join(data_dir, 'glg_tte_n9_bn090131090_v00.fit')
    trigtime = 255060563.149072
    numchans = 128
    erange = (4.389729, 2000.0)
    gti = (255060537.657974, 255060863.884596)
    gti_rel = (gti[0]-trigtime, gti[1]-trigtime)
    trange = (255060537.657974, 255060863.884596)
    trange_rel = (trange[0]-trigtime, trange[1]-trigtime)

    def test_attributes(self):
        tte = TTE.open(self.filename)
        trigtime = tte.trigtime
        self.assertAlmostEqual(trigtime, self.trigtime, places=6)
        self.assertEqual(tte.is_gbm_file, True)
        self.assertEqual(tte.id, '090131090')
        self.assertEqual(tte.filename, os.path.basename(self.filename))
        self.assertEqual(tte.is_trigger, True)
        self.assertEqual(tte.detector, 'n9')
        self.assertEqual(tte.datatype, 'TTE')
        self.assertAlmostEqual(tte.energy_range[0], self.erange[0], places=4)
        self.assertAlmostEqual(tte.energy_range[1], self.erange[1], places=4)
        self.assertEqual(len(tte.gti), 1)
        self.assertAlmostEqual(tte.gti[0][0], self.gti_rel[0], places=6)
        self.assertAlmostEqual(tte.gti[0][1], self.gti_rel[1], places=6)
        self.assertCountEqual(tte.headers.keys(), ['PRIMARY', 'EBOUNDS', 
                                                     'EVENTS', 'GTI'])
        self.assertEqual(tte.numchans, self.numchans)
        self.assertAlmostEqual(tte.time_range[0], self.trange_rel[0], places=6)
        self.assertAlmostEqual(tte.time_range[1], self.trange_rel[1], places=6)

    def test_get_exposure(self):
        tte = TTE.open(self.filename)
        exposure = tte.get_exposure(time_ranges=(0.0, 10.0))
        self.assertAlmostEqual(exposure, 10.0, delta=0.1)        
        exposure = tte.get_exposure(time_ranges=[(0.0, 10.0)])
        self.assertAlmostEqual(exposure, 10.00, delta=0.1)        
        exposure = tte.get_exposure(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
        self.assertAlmostEqual(exposure, 20.0, delta=0.2)        

    def test_to_phaii(self):
        tte = TTE.open(self.filename)
        pha1 = tte.to_phaii(bin_by_time, 1.024)
        self.assertEqual(pha1.time_range[0], tte.time_range[0])
        self.assertCountEqual(pha1.energy_range, tte.energy_range)
        self.assertEqual(pha1.data.numtimes, 319)

        pha2 = tte.to_phaii(bin_by_time, 1.024, time_range=(-10.0, 100.0))
        self.assertEqual(pha2.time_range[0], -10.0)
        self.assertCountEqual(pha2.energy_range, tte.energy_range)
        self.assertEqual(pha2.data.numtimes, 108)

        pha3 = tte.to_phaii(bin_by_time, 1.024, energy_range=(50.0, 300.0))
        self.assertAlmostEqual(pha3.time_range[0], tte.time_range[0], places=1)
        self.assertCountEqual(pha3.energy_range, tte.energy_range)
        self.assertEqual(pha3.data.numtimes, 319)

        pha4 = tte.to_phaii(bin_by_time, 1.024, channel_range=(50, 80))
        self.assertAlmostEqual(pha4.time_range[0], tte.time_range[0], places=1)
        self.assertCountEqual(pha4.energy_range, tte.energy_range)
        self.assertEqual(pha4.data.numtimes, 319)

    def test_to_pha(self):
        tte = TTE.open(self.filename)
        pha1 = tte.to_pha()
        self.assertCountEqual(pha1.time_range, tte.time_range)
        self.assertCountEqual(pha1.energy_range, tte.energy_range)
        self.assertCountEqual(pha1.valid_channels, np.arange(self.numchans))
        
        pha2 = tte.to_pha(energy_range=(50.0, 300.0))
        self.assertCountEqual(pha2.time_range, tte.time_range)
        self.assertCountEqual(pha2.energy_range, tte.energy_range)
        self.assertCountEqual(pha2.valid_channels, np.arange(33,85))

        pha3 = tte.to_pha(channel_range=(50, 80))
        self.assertCountEqual(pha3.time_range, tte.time_range)
        self.assertCountEqual(pha3.energy_range, tte.energy_range)
        self.assertCountEqual(pha3.valid_channels, np.arange(50,81))

        pha4 = tte.to_pha(time_ranges=[(-10.0,10.0), (20.0, 30.0)])
        self.assertAlmostEqual(pha4.time_range[0], -10.0, places=1)
        self.assertAlmostEqual(pha4.time_range[1], 30.0, places=1)
        self.assertCountEqual(pha4.energy_range, tte.energy_range)
        self.assertCountEqual(pha4.valid_channels, np.arange(self.numchans))

    def test_write(self):
        tte = TTE.open(self.filename)
        new_file = os.path.join(data_dir, 'glg_tte_n9_bn090131090_v01.fit')
        tte.write(os.path.dirname(self.filename), 
                  filename='glg_tte_n9_bn090131090_v01.fit')    
        tte2 = TTE.open(new_file)
        tte2.data.sort('TIME')
        os.remove(new_file)

        self.assertEqual(tte2.data.size, tte.data.size)
        self.assertCountEqual(tte2.data.pha, tte.data.pha)
        self.assertCountEqual(tte2.data.emin, tte.data.emin)
        self.assertCountEqual(tte2.data.emax, tte.data.emax)
        self.assertCountEqual(tte2.data.time, tte.data.time)

    def test_slice_time(self):
        tte = TTE.open(self.filename)
        tte2 = tte.slice_time((-10.0, 10.0))
        self.assertEqual(tte2.gti, [tte2.time_range])
        
        tte3 = tte.slice_time([(-10.0, 10.0), (20.0, 30.0)])        
        self.assertEqual(len(tte3.gti), 2)

    def test_slice_energy(self):
        tte = TTE.open(self.filename)
        tte2 = tte.slice_energy((50.0, 250.0))
        self.assertAlmostEqual(tte2.gti[0][0], tte2.time_range[0], delta=0.1)
        self.assertAlmostEqual(tte2.gti[0][1], tte2.time_range[1], delta=0.1)
        
        tte3 = tte.slice_energy([(50.0, 250.0), (500.0, 2000.0)])        
        self.assertAlmostEqual(tte2.gti[0][0], tte2.time_range[0], delta=0.1)
        self.assertAlmostEqual(tte2.gti[0][1], tte2.time_range[1], delta=0.1)
    
    def test_merge(self):
        tte = TTE.open(self.filename)
        tte2 = tte.slice_time((-10.0, 10.0))
        tte3 = TTE.merge((tte, tte2))
        
        self.assertEqual(tte.time_range[0], tte3.time_range[0])
        self.assertEqual(tte.time_range[1], tte3.time_range[1])
        self.assertEqual(tte.data.size, tte3.data.size)
    
    def test_rebin_energy(self):
        tte = TTE.open(self.filename)
        tte2 = tte.rebin_energy(combine_by_factor, 4)
        
        self.assertEqual(tte.time_range[0], tte2.time_range[0])
        self.assertEqual(tte.time_range[1], tte2.time_range[1])
        self.assertEqual(tte.data.size, tte2.data.size)
        self.assertEqual(tte.data.channel_range[0], 0)
        self.assertEqual(tte.data.channel_range[1], 31)
       
    def test_errors(self):
        with  self.assertRaises(IOError):
            tte = TTE.open('stephenking.fit')
        
        tte = TTE.open(self.filename)
        with self.assertRaises(AssertionError):
            tte.to_pha(energy_range=(30.0, 10.0))
            tte.to_pha(channel_range=(30, 10))
            tte.to_pha(time_ranges=(100.0, 10.0))

