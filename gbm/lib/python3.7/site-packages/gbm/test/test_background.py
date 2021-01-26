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
import os
import numpy as np
from unittest import TestCase

from gbm.background.binned import Polynomial
from gbm.background.unbinned import NaivePoisson
from gbm.background.background import BackgroundFitter, BackgroundRates, BackgroundSpectrum
from gbm.data.phaii import Cspec, TTE

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestPolynomialBackground(TestCase):
    counts = np.array([[78.0, 58.0, 40.0, 26.0, 14.0, 6.0, 2.0, 0.0, 2.0, 6.0],
                       [6.0, 2.0, 0.0, 2.0, 6.0, 14.0, 26.0, 40.0, 58.0, 78.0]]).T
    exposure = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
    edges = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    tstart = edges[:-1]
    tstop = edges[1:]
    
    vals = np.array([[39.0, 29.0, 20.0, 13.0, 7.0, 3.0, 1.0, 0.0, 1.0, 3.0], 
                     [3.0, 1.0, 0.0, 1.0, 3.0, 7.0, 13.0, 20.0, 29.0, 39.0]]).T
    
    def test_fit(self): 
        bkgd = Polynomial(self.counts, self.tstart, self.tstop, self.exposure)
        bkgd.fit(order=2)
        self.assertEqual(bkgd.statistic_name, 'chisq')
        self.assertCountEqual(bkgd.dof, np.array([6, 6]))
        self.assertAlmostEqual(bkgd.statistic[0], 0.23, delta=0.01)
        self.assertAlmostEqual(bkgd.statistic[1], 0.23, delta=0.01)
    
    def test_interpolate(self):
        bkgd = Polynomial(self.counts, self.tstart, self.tstop, self.exposure)
        bkgd.fit(order=2)
        rates, rate_uncert = bkgd.interpolate(self.tstart, self.tstop)
        for i in range(10):
            self.assertAlmostEqual(rates[i,0], self.vals[i,0], delta=0.5)
            self.assertAlmostEqual(rates[i,1], self.vals[i,1], delta=0.5)
    
class TestNaivePoissonBackground(TestCase):
    times = [np.array([0., 1.14, 1.22, 1.28, 1.76, 3.45, 4.29, 4.78, 4.75, 5.42,
                      5.97, 7.40, 7.61, 7.98, 8.10, 8.16, 10.18, 10.13, 
                      13.22, 14.03])]*2
    
    def test_fit(self):
        bkgd = NaivePoisson(self.times)
        bkgd.fit(window_width=5.0, fast=True)
        bkgd.fit(window_width=5.0, fast=False)
    
    def test_interpolate(self):
        bkgd = NaivePoisson(self.times)
        bkgd.fit(window_width=5.0, fast=True)
        x = np.linspace(0.0, 14.0, 15)
        rates, uncert = bkgd.interpolate(x[:-1], x[1:])
        for i in range(14):
            self.assertAlmostEqual(rates[i,0], 1.0, delta=2.)
            self.assertAlmostEqual(rates[i,1], 1.0, delta=2.)

        bkgd.fit(window_width=5.0, fast=False)
        x = np.linspace(0.0, 14.0, 15)
        rates, uncert = bkgd.interpolate(x[:-1], x[1:])
        for i in range(14):
            self.assertAlmostEqual(rates[i,0], 1.0, delta=2.)
            self.assertAlmostEqual(rates[i,1], 1.0, delta=2.)

class TestBackgroundFitterBinned(TestCase):
    file = os.path.join(data_dir, 'glg_cspec_b0_bn120415958_v00.pha')
    
    def test_attributes(self):
        cspec = Cspec.open(self.file)
        b = BackgroundFitter.from_phaii(cspec, Polynomial, 
                                        time_ranges=[(-100, -10.0), (100.0, 200.0)])
        b.fit(order=3)
        self.assertEqual(b.method, 'Polynomial')
        self.assertEqual(b.type, 'binned')
        self.assertEqual(b.statistic_name, 'chisq')
        self.assertEqual(b.dof.size, 128)
        self.assertEqual(b.statistic.size, 128)
        self.assertAlmostEqual(b.livetime, 190.0, delta=5.0)
        self.assertEqual(b.parameters, {'order': 3})
    
    def test_interpolate_bins(self):
        cspec = Cspec.open(self.file)
        b = BackgroundFitter.from_phaii(cspec, Polynomial, 
                                        time_ranges=[(-100, -10.0), (100.0, 200.0)])
        b.fit(order=2)
        x = np.linspace(-10.0, 11.0, 21)
        tstart = x[:-1]
        tstop = x[1:]
        brates = b.interpolate_bins(tstart, tstop)
        self.assertEqual(brates.__class__.__name__, 'BackgroundRates')

class TestBackgroundFitterUnbinned(TestCase):
    file = os.path.join(data_dir, 'glg_tte_n9_bn090131090_v00.fit')

    def test_attributes(self):
        tte = TTE.open(self.file)
        b = BackgroundFitter.from_tte(tte, NaivePoisson, 
                                        time_ranges=[(-20, 200.0)])
        b.fit(window_width=30.0, fast=True)
        self.assertEqual(b.method, 'NaivePoisson')
        self.assertEqual(b.type, 'unbinned')
        self.assertAlmostEqual(b.livetime, 220.0, delta=1.0)
        self.assertEqual(b.parameters, {'fast': True, 'window_width': 30.0})

    def test_attributes(self):
        tte = TTE.open(self.file)
        b = BackgroundFitter.from_tte(tte, NaivePoisson, 
                                        time_ranges=[(-20, 200.0)])
        b.fit(window_width=30.0, fast=True)
        x = np.linspace(-20.0, 200.0, 221)
        tstart = x[:-1]
        tstop = x[1:]
        brates = b.interpolate_bins(tstart, tstop)
        self.assertEqual(brates.__class__.__name__, 'BackgroundRates')
     
class TestBackgroundRates(TestCase):
    file = os.path.join(data_dir, 'glg_cspec_b0_bn120415958_v00.pha')
    
    def test_attributes(self):
        cspec = Cspec.open(self.file)
        b = BackgroundFitter.from_phaii(cspec, Polynomial, 
                                        time_ranges=[(-100, -10.0), (100.0, 200.0)])
        b.fit(order=3)
        
        b.fit(order=2)
        tstart = cspec.data.tstart
        mask = (tstart >= -10.0) & (tstart < 11.0)
        tstart = tstart[mask]
        tstop = cspec.data.tstop[mask]
        brates = b.interpolate_bins(tstart, tstop)

        self.assertEqual(brates.numchans, 128)
        self.assertEqual(brates.numtimes, 12)
        self.assertEqual(brates.size, (12, 128))
        self.assertCountEqual(brates.emin, cspec.data.emin)
        self.assertCountEqual(brates.emax, cspec.data.emax)
        self.assertAlmostEqual(np.sum(brates.exposure), 17.0, delta=1.0)
        
        bspec = brates.integrate_time()
        self.assertEqual(bspec.__class__.__name__, 'BackgroundSpectrum')
        bak = brates.to_bak()
        self.assertEqual(bak.__class__.__name__, 'BAK')
        
        brates2 = brates.integrate_energy()
        self.assertEqual(brates2.counts.size, brates2.rates.size)
        self.assertEqual(brates2.count_uncertainty.size, brates2.rates.size)
        

class TestBackgroundSpectrum(TestCase):
    file = os.path.join(data_dir, 'glg_cspec_b0_bn120415958_v00.pha')
    
    def test_attributes(self):
        cspec = Cspec.open(self.file)
        b = BackgroundFitter.from_phaii(cspec, Polynomial, 
                                        time_ranges=[(-100, -10.0), (100.0, 200.0)])
        b.fit(order=3)
        
        b.fit(order=2)
        x = np.linspace(-10.0, 11.0, 21)
        tstart = x[:-1]
        tstop = x[1:]
        brates = b.interpolate_bins(tstart, tstop)
        bspec = brates.integrate_time(tstart=-5.0, tstop=5.0)
        self.assertEqual(bspec.size, 128)
        self.assertCountEqual(bspec.lo_edges, cspec.data.emin)
        self.assertCountEqual(bspec.hi_edges, cspec.data.emax)
    
    def test_slice(self):
        cspec = Cspec.open(self.file)
        b = BackgroundFitter.from_phaii(cspec, Polynomial, 
                                        time_ranges=[(-100, -10.0), (100.0, 200.0)])
        b.fit(order=3)
        x = np.linspace(-10.0, 11.0, 21)
        tstart = x[:-1]
        tstop = x[1:]
        brates = b.interpolate_bins(tstart, tstop)
        bspec = brates.integrate_time(tstart=-5.0, tstop=5.0)
        bspec_sliced = bspec.slice(100.0, 300.0)
        self.assertCountEqual(bspec.counts[:4], bspec_sliced.counts)
        self.assertCountEqual(bspec.rates[:4], bspec_sliced.rates)
        self.assertCountEqual(bspec.exposure[:4], bspec_sliced.exposure)
        self.assertCountEqual(bspec.lo_edges[:4], bspec_sliced.lo_edges)
        self.assertCountEqual(bspec.hi_edges[:4], bspec_sliced.hi_edges)
        self.assertCountEqual(bspec.rate_uncertainty[:4], bspec_sliced.rate_uncertainty)
              
if __name__ == '__main__':
    unittest.main()
      
        
        