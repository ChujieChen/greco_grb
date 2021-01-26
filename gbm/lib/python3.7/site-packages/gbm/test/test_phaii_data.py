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

from gbm.data.phaii import Ctime, Cspec
from gbm.binning.binned import combine_by_factor

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestCtime(TestCase):
    filename = os.path.join(data_dir, 'glg_ctime_nb_bn120415958_v00.pha')
    trigtime = 356223561.133346
    numchans = 8
    erange = (4.323754, 2000.0)
    gti = (356222661.790904, 356224561.991216)
    gti_rel = (gti[0]-trigtime, gti[1]-trigtime)
    trange = (356222661.790904, 356224561.991216)
    trange_rel = (trange[0]-trigtime, trange[1]-trigtime)
    
    def test_attributes(self):
        ctime = Ctime.open(self.filename)
        trigtime = ctime.trigtime
        self.assertAlmostEqual(trigtime, self.trigtime, places=6)
        self.assertEqual(ctime.is_gbm_file, True)
        self.assertEqual(ctime.id, '120415958')
        self.assertEqual(ctime.filename, os.path.basename(self.filename))
        self.assertEqual(ctime.is_trigger, True)
        self.assertEqual(ctime.detector, 'nb')
        self.assertEqual(ctime.datatype, 'CTIME')
        self.assertAlmostEqual(ctime.energy_range[0], self.erange[0], places=6)
        self.assertAlmostEqual(ctime.energy_range[1], self.erange[1], places=6)
        self.assertEqual(len(ctime.gti), 1)
        self.assertAlmostEqual(ctime.gti[0][0], self.gti_rel[0], places=6)
        self.assertAlmostEqual(ctime.gti[0][1], self.gti_rel[1], places=6)
        self.assertCountEqual(ctime.headers.keys(), ['PRIMARY', 'EBOUNDS', 
                                                     'SPECTRUM', 'GTI'])
        self.assertEqual(ctime.numchans, self.numchans)
        self.assertAlmostEqual(ctime.time_range[0], self.trange_rel[0], places=6)
        self.assertAlmostEqual(ctime.time_range[1], self.trange_rel[1], places=6)
        self.assertEqual(len(ctime.data.contiguous_time_bins()), 1)
        self.assertEqual(len(ctime.data.contiguous_energy_bins()), 1)
              
    def test_get_exposure(self):
        ctime = Ctime.open(self.filename)
        exposure = ctime.get_exposure(time_ranges=(0.0, 10.0))
        self.assertAlmostEqual(exposure, 10.01, delta=0.1)        
        exposure = ctime.get_exposure(time_ranges=[(0.0, 10.0)])
        self.assertAlmostEqual(exposure, 10.01, delta=0.1)        
        exposure = ctime.get_exposure(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
        self.assertAlmostEqual(exposure, 20.02, delta=0.1)        

    def test_to_lightcurve(self):
        ctime = Ctime.open(self.filename)
        lc1 = ctime.to_lightcurve()       
        self.assertAlmostEqual(lc1.range[0], self.trange[0]-self.trigtime, delta=0.1)
        self.assertAlmostEqual(lc1.range[1], self.trange[1]-self.trigtime, delta=0.1)
        self.assertEqual(lc1.size, 14433)
        
        lc2 = ctime.to_lightcurve(time_range=(-100.0, 100.0))
        self.assertAlmostEqual(lc2.range[0], 356223461.035746-self.trigtime, delta=0.1)
        self.assertAlmostEqual(lc2.range[1], 356223661.16687-self.trigtime, delta=0.1)
        self.assertEqual(lc2.size, 1955)

        lc3 = ctime.to_lightcurve(energy_range=(50.0, 300.0))       
        self.assertAlmostEqual(lc3.range[0], self.trange[0]-self.trigtime, delta=0.1)
        self.assertAlmostEqual(lc3.range[1], self.trange[1]-self.trigtime, delta=0.1)
        self.assertEqual(lc3.size, 14433)
        
        lc4 = ctime.to_lightcurve(channel_range=(1, 6))       
        self.assertAlmostEqual(lc4.range[0], self.trange[0]-self.trigtime, delta=0.1)
        self.assertAlmostEqual(lc4.range[1], self.trange[1]-self.trigtime, delta=0.1)
        self.assertEqual(lc4.size, 14433)

    def test_to_spectrum(self):
        ctime = Ctime.open(self.filename)
        lc1 = ctime.to_spectrum()       
        self.assertAlmostEqual(lc1.range[0], self.erange[0], delta=0.1)
        self.assertAlmostEqual(lc1.range[1], self.erange[1], delta=0.1)
        self.assertEqual(lc1.size, 8)
        
        lc2 = ctime.to_spectrum(time_range=(-100.0, 100.0))
        self.assertAlmostEqual(lc2.range[0], self.erange[0], delta=0.1)
        self.assertAlmostEqual(lc2.range[1], self.erange[1], delta=0.1)
        self.assertEqual(lc2.size, 8)

        lc3 = ctime.to_spectrum(energy_range=(50.0, 300.0))       
        self.assertAlmostEqual(lc3.range[0], 49.60019, delta=0.1)
        self.assertAlmostEqual(lc3.range[1], 538.144, delta=0.1)
        self.assertEqual(lc3.size, 3)
        
        lc4 = ctime.to_spectrum(channel_range=(3, 4))       
        self.assertAlmostEqual(lc4.range[0], 49.60019, delta=0.1)
        self.assertAlmostEqual(lc4.range[1], 290.4606, delta=0.1)
        self.assertEqual(lc4.size, 2)
    
    def test_to_pha(self):
        ctime = Ctime.open(self.filename)
        pha1 = ctime.to_pha()
        self.assertCountEqual(pha1.time_range, ctime.time_range)      
        self.assertCountEqual(pha1.gti, ctime.gti)
        self.assertEqual(pha1.trigtime, ctime.trigtime)
        self.assertAlmostEqual(pha1.energy_range[0], self.erange[0], delta=0.1)
        self.assertAlmostEqual(pha1.energy_range[1], self.erange[1], delta=0.1)
        self.assertEqual(pha1.numchans, self.numchans)
        self.assertCountEqual(pha1.valid_channels, np.arange(self.numchans))
        self.assertAlmostEqual(pha1.exposure, np.sum(ctime.data.exposure), delta=0.1)
        
        pha2 = ctime.to_pha(channel_range=(3,4))       
        self.assertAlmostEqual(pha2.energy_range[0], self.erange[0], delta=0.1)
        self.assertAlmostEqual(pha2.energy_range[1], self.erange[1], delta=0.1)
        self.assertEqual(pha2.numchans, self.numchans)
        self.assertCountEqual(pha2.valid_channels, np.array([3,4]))
        
        pha3 = ctime.to_pha(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
        self.assertAlmostEqual(pha3.time_range[0], 0.0, delta=0.1)   
        self.assertAlmostEqual(pha3.time_range[1], 30.0, delta=0.1)   
        self.assertAlmostEqual(pha3.exposure, 20., delta=0.1)
  
    def test_write(self):
        ctime = Ctime.open(self.filename)
        new_file = os.path.join(data_dir, 'glg_ctime_nb_bn120415958_v01.pha')
        ctime.write(os.path.dirname(self.filename), 
                    filename='glg_ctime_nb_bn120415958_v01.pha')    
        ctime2 = Ctime.open(new_file)
        os.remove(new_file)
        
        self.assertCountEqual(ctime2.data.tstart, ctime.data.tstart)
        self.assertCountEqual(ctime2.data.tstop, ctime.data.tstop)
        self.assertCountEqual(ctime2.data.emin, ctime.data.emin)
        self.assertCountEqual(ctime2.data.emax, ctime.data.emax)
        self.assertCountEqual(ctime2.data.exposure, ctime.data.exposure)
        for i in range(ctime.numchans):
            self.assertCountEqual(ctime2.data.counts[:,i],
                                  ctime.data.counts[:,i])
    
    def test_slice_time(self):
        ctime = Ctime.open(self.filename)
        ctime2 = ctime.slice_time((-10.0, 10.0))
        self.assertEqual(ctime2.gti, [ctime2.time_range])
        self.assertEqual(ctime2.data.numtimes, 198)
        
        ctime3 = ctime.slice_time([(-10.0, 10.0), (20.0, 30.0)])        
        self.assertEqual(len(ctime3.gti), 2)
        self.assertEqual(ctime3.data.numtimes, 355)
        
    def test_slice_energy(self):
        ctime = Ctime.open(self.filename)
        ctime2 = ctime.slice_energy((50.0, 250.0))
        self.assertAlmostEqual(ctime2.energy_range[0], 50.0, delta=1.0)
        self.assertAlmostEqual(ctime2.energy_range[1], 250.0, delta=50.0)
        self.assertEqual(ctime2.data.numchans, 2)
        
        ctime3 = ctime.slice_energy([(50.0, 250.0), (500.0, 2000.0)])        
        self.assertAlmostEqual(ctime3.energy_range[0], 50.0, delta=1.0)
        self.assertAlmostEqual(ctime3.energy_range[1], 2000.0, delta=50.0)
        self.assertEqual(ctime3.data.numchans, 5)
        
    def test_rebin_time(self):
        ctime = Ctime.open(self.filename)
        ctime2 = ctime.rebin_time(combine_by_factor, 4)
        self.assertAlmostEqual(ctime2.gti[0][0], ctime.gti[0][0], delta=1.0)
        self.assertAlmostEqual(ctime2.gti[0][1], ctime.gti[0][1], delta=1.0)
        self.assertEqual(ctime2.data.numtimes, 3608)
        
        ctime3 = ctime.rebin_time(combine_by_factor, 4, time_range=(-100.0, 100.0))       
        self.assertEqual(ctime3.gti, ctime.gti)        
        self.assertEqual(ctime3.data.numtimes, 12968)

    def test_rebin_energy(self):
        ctime = Ctime.open(self.filename)
        ctime2 = ctime.rebin_energy(combine_by_factor, 2)
        self.assertEqual(ctime2.energy_range, ctime.energy_range)
        self.assertEqual(ctime2.data.numchans, 4)
        
        ctime3 = ctime.rebin_energy(combine_by_factor, 4, energy_range=(30.0, 300.0))  
        self.assertEqual(ctime3.energy_range, ctime.energy_range)
        self.assertEqual(ctime3.data.numchans, 6)
    
    def test_errors(self):
        with  self.assertRaises(IOError):
            ctime = Ctime.open('wakawaka.pha')
        
        ctime = Ctime.open(self.filename)
        with self.assertRaises(AssertionError):
            ctime.to_lightcurve(time_range=(100.0, 10.0))
            ctime.to_lightcurve(energy_range=(30.0, 10.0))
            ctime.to_lightcurve(channel_range=(30, 10))
            ctime.to_spectrum(time_range=(100.0, 10.0))
            ctime.to_spectrum(energy_range=(30.0, 10.0))
            ctime.to_spectrum(channel_range=(30, 10))
        

class TestCspec(TestCase):
    filename = os.path.join(data_dir, 'glg_cspec_b0_bn120415958_v00.pha')
    trigtime = 356223561.133346
    numchans = 128
    erange = (114.5423, 50000.0)
    gti = [(356219559.260874, 356221211.030666), (356222657.950904, 356227367.05709)]
    gti_rel = [(gti[0][0]-trigtime, gti[0][1]-trigtime), 
               (gti[1][0]-trigtime, gti[1][1]-trigtime)]
    trange = (356219559.260874, 356227367.057090 )
    trange_rel = (trange[0]-trigtime, trange[1]-trigtime)

    def test_attributes(self):
        cspec = Cspec.open(self.filename)
        trigtime = cspec.trigtime
        self.assertAlmostEqual(trigtime, self.trigtime, places=6)
        self.assertEqual(cspec.is_gbm_file, True)
        self.assertEqual(cspec.id, '120415958')
        self.assertEqual(cspec.filename, os.path.basename(self.filename))
        self.assertEqual(cspec.is_trigger, True)
        self.assertEqual(cspec.detector, 'b0')
        self.assertEqual(cspec.datatype, 'CSPEC')
        self.assertAlmostEqual(cspec.energy_range[0], self.erange[0], places=4)
        self.assertAlmostEqual(cspec.energy_range[1], self.erange[1], places=4)
        self.assertEqual(len(cspec.gti), 2)
        self.assertAlmostEqual(cspec.gti[0][0], self.gti_rel[0][0], places=6)
        self.assertAlmostEqual(cspec.gti[0][1], self.gti_rel[0][1], places=6)
        self.assertAlmostEqual(cspec.gti[1][0], self.gti_rel[1][0], places=6)
        self.assertAlmostEqual(cspec.gti[1][1], self.gti_rel[1][1], places=6)
        self.assertCountEqual(cspec.headers.keys(), ['PRIMARY', 'EBOUNDS', 
                                                     'SPECTRUM', 'GTI'])
        self.assertEqual(cspec.numchans, self.numchans)
        self.assertAlmostEqual(cspec.time_range[0], self.trange_rel[0], places=6)
        self.assertAlmostEqual(cspec.time_range[1], self.trange_rel[1], places=6)
        self.assertEqual(len(cspec.data.contiguous_time_bins()), 3)
        self.assertEqual(len(cspec.data.contiguous_energy_bins()), 1)

    def test_to_lightcurve(self):
        cspec = Cspec.open(self.filename)
        lc1 = cspec.to_lightcurve()       
        self.assertAlmostEqual(lc1.range[0], self.trange[0]-self.trigtime, places=6)
        self.assertAlmostEqual(lc1.range[1], self.trange[1]-self.trigtime, places=6)
        self.assertEqual(lc1.size, 1994)
        
        lc2 = cspec.to_lightcurve(time_range=(-100.0, 100.0))
        self.assertAlmostEqual(lc2.range[0], 356223460.77972996-self.trigtime, places=6)
        self.assertAlmostEqual(lc2.range[1], 356223661.48687-self.trigtime, places=6)
        self.assertEqual(lc2.size, 122)

        lc3 = cspec.to_lightcurve(energy_range=(500.0, 3000.0))       
        self.assertAlmostEqual(lc3.range[0], self.trange[0]-self.trigtime, places=6)
        self.assertAlmostEqual(lc3.range[1], self.trange[1]-self.trigtime, places=6)
        self.assertEqual(lc3.size, 1994)
        
        lc4 = cspec.to_lightcurve(channel_range=(50, 80))       
        self.assertAlmostEqual(lc4.range[0], self.trange[0]-self.trigtime, places=6)
        self.assertAlmostEqual(lc4.range[1], self.trange[1]-self.trigtime, places=6)
        self.assertEqual(lc4.size, 1994)

    def test_to_spectrum(self):
        cspec = Cspec.open(self.filename)
        lc1 = cspec.to_spectrum()       
        self.assertAlmostEqual(lc1.range[0], self.erange[0], places=4)
        self.assertAlmostEqual(lc1.range[1], self.erange[1], places=4)
        self.assertEqual(lc1.size, self.numchans)
        
        lc2 = cspec.to_spectrum(time_range=(-100.0, 100.0))
        self.assertAlmostEqual(lc2.range[0], self.erange[0], places=4)
        self.assertAlmostEqual(lc2.range[1], self.erange[1], places=4)
        self.assertEqual(lc2.size, self.numchans)

        lc3 = cspec.to_spectrum(energy_range=(500.0, 3000.0))       
        self.assertAlmostEqual(lc3.range[0], 488.1802, places=3)
        self.assertAlmostEqual(lc3.range[1], 3087.304, places=3)
        self.assertEqual(lc3.size, 43)
        
        lc4 = cspec.to_spectrum(channel_range=(50, 80))       
        self.assertAlmostEqual(lc4.range[0], 2894.153, places=3)
        self.assertAlmostEqual(lc4.range[1], 7464.066, places=3)
        self.assertEqual(lc4.size, 31)

    def test_to_pha(self):
        cspec = Cspec.open(self.filename)
        pha1 = cspec.to_pha()
        self.assertCountEqual(pha1.time_range, cspec.time_range)      
        self.assertCountEqual(pha1.gti[0], cspec.time_range)
        self.assertEqual(pha1.trigtime, cspec.trigtime)
        self.assertAlmostEqual(pha1.energy_range[0], self.erange[0], places=4)
        self.assertAlmostEqual(pha1.energy_range[1], self.erange[1], places=4)
        self.assertEqual(pha1.numchans, self.numchans)
        self.assertCountEqual(pha1.valid_channels, np.arange(self.numchans))
        self.assertAlmostEqual(pha1.exposure, np.sum(cspec.data.exposure), places=5)
        
        pha2 = cspec.to_pha(channel_range=(50,80))       
        self.assertAlmostEqual(pha2.energy_range[0], self.erange[0], places=4)
        self.assertAlmostEqual(pha2.energy_range[1], self.erange[1], places=4)
        self.assertEqual(pha2.numchans, self.numchans)
        self.assertCountEqual(pha2.valid_channels, np.arange(50, 81))
        
        pha3 = cspec.to_pha(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
        self.assertAlmostEqual(pha3.time_range[0], -2.048, places=1)   
        self.assertAlmostEqual(pha3.time_range[1], 31.0, places=0)   
        self.assertAlmostEqual(pha3.exposure, 23.3, places=0)

    def test_write(self):
        cspec = Cspec.open(self.filename)
        new_file = os.path.join(data_dir, 'glg_cspec_b0_bn120415958_v01.pha')
        cspec.write(os.path.dirname(self.filename), 
                    filename='glg_cspec_b0_bn120415958_v01.pha')    
        cspec2 = Cspec.open(new_file)
        os.remove(new_file)
        
        self.assertCountEqual(cspec2.data.tstart, cspec.data.tstart)
        self.assertCountEqual(cspec2.data.tstop, cspec.data.tstop)
        self.assertCountEqual(cspec2.data.emin, cspec.data.emin)
        self.assertCountEqual(cspec2.data.emax, cspec.data.emax)
        self.assertCountEqual(cspec2.data.exposure, cspec.data.exposure)
        for i in range(cspec.numchans):
            self.assertCountEqual(cspec2.data.counts[:,i],
                                  cspec.data.counts[:,i])

    def test_slice_time(self):
        cspec = Cspec.open(self.filename)
        cspec2 = cspec.slice_time((-10.0, 10.0))
        self.assertEqual(cspec2.gti, [cspec2.time_range])
        self.assertEqual(cspec2.data.numtimes, 12)
        
        cspec3 = cspec.slice_time([(-10.0, 10.0), (20.0, 30.0)])        
        self.assertEqual(len(cspec3.gti), 2)
        self.assertEqual(cspec3.data.numtimes, 23)
        
    def test_slice_energy(self):
        cspec = Cspec.open(self.filename)
        cspec2 = cspec.slice_energy((300.0, 1000.0))
        self.assertAlmostEqual(cspec2.energy_range[0], 300.0, delta=50.0)
        self.assertAlmostEqual(cspec2.energy_range[1], 1000.0, delta=50.0)
        self.assertEqual(cspec2.data.numchans, 19)
        
        cspec3 = cspec.slice_energy([(300.0, 1000.0), (5000.0, 10000.0)])        
        self.assertAlmostEqual(cspec3.energy_range[0], 300.0, delta=50.0)
        self.assertAlmostEqual(cspec3.energy_range[1], 10000.0, delta=200.0)
        self.assertEqual(cspec3.data.numchans, 43)
        
    def test_rebin_time(self):
        cspec = Cspec.open(self.filename)
        cspec2 = cspec.rebin_time(combine_by_factor, 4)
        self.assertAlmostEqual(cspec2.gti[0][0], cspec.gti[0][0], delta=1.0)
        self.assertAlmostEqual(cspec2.gti[-1][1], cspec.gti[-1][1], delta=10.0)
        self.assertEqual(cspec2.data.numtimes, 497)
        
        cspec3 = cspec.rebin_time(combine_by_factor, 4, time_range=(-100.0, 100.0))       
        self.assertEqual(cspec3.data.numtimes, 1904)

    def test_rebin_energy(self):
        cspec = Cspec.open(self.filename)
        cspec2 = cspec.rebin_energy(combine_by_factor, 4)
        self.assertEqual(cspec2.energy_range, cspec.energy_range)
        self.assertEqual(cspec2.data.numchans, 32)
        
        cspec3 = cspec.rebin_energy(combine_by_factor, 4, energy_range=(500., 2000.0))       
        self.assertEqual(cspec3.energy_range, cspec.energy_range)
        self.assertEqual(cspec3.data.numchans, 107)

    def test_errors(self):
        with  self.assertRaises(IOError):
            cspec = Cspec.open('michaelscott.pha')
        
        cspec = Ctime.open(self.filename)
        with self.assertRaises(AssertionError):
            cspec.to_lightcurve(time_range=(100.0, 10.0))
            cspec.to_lightcurve(energy_range=(30.0, 10.0))
            cspec.to_lightcurve(channel_range=(30, 10))
            cspec.to_spectrum(time_range=(100.0, 10.0))
            cspec.to_spectrum(energy_range=(30.0, 10.0))
            cspec.to_spectrum(channel_range=(30, 10))

if __name__ == '__main__':
    unittest.main()