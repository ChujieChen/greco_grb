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

from gbm.simulate.generators import *
from gbm.simulate.profiles import *
from gbm.simulate import PhaSimulator, TteSourceSimulator, TteBackgroundSimulator
from gbm.data import BAK, RSP, PHA, PHAII, TTE
from gbm.spectra.functions import Band

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestGenerators(TestCase):
    bkgd = BAK.open(os.path.join(data_dir, 'glg_tte_n0_bn160509374_xspec_v00.bak'))
    bkgd = bkgd.data
    rsp = RSP.open(os.path.join(data_dir, 'glg_cspec_n4_bn120415958_v00.rsp2'))
    fxn = Band()
    params = (1.0, 300.0, -1.0, -2.5)
    exposure = 2.048
    
    def test_poisson_background(self): 
        gen = PoissonBackgroundGenerator(self.bkgd)
        dev = next(gen)
        self.assertEqual(dev.size, self.bkgd.size)
        self.assertCountEqual(dev.lo_edges, self.bkgd.lo_edges)
        self.assertCountEqual(dev.hi_edges, self.bkgd.hi_edges)
        self.assertCountEqual(dev.exposure, self.bkgd.exposure)
        devs = [next(gen) for i in range(10)]

    def test_variable_poisson_background(self): 
        gen = VariablePoissonBackground(self.bkgd)
        dev = next(gen)
        self.assertEqual(dev.size, self.bkgd.size)
        self.assertCountEqual(dev.lo_edges, self.bkgd.lo_edges)
        self.assertCountEqual(dev.hi_edges, self.bkgd.hi_edges)
        self.assertCountEqual(dev.exposure, self.bkgd.exposure)
        devs = [next(gen) for i in range(10)]
        
        gen.amp = 10.0
        dev = next(gen)
        ratio = dev.counts/devs[0].counts
        for i in range(self.bkgd.size):
            if not np.isnan(ratio[i]):
                self.assertTrue(ratio[i] > 1.0)

        gen.amp = 0.1
        dev = next(gen)
        ratio = dev.counts/devs[0].counts
        for i in range(self.bkgd.size):
            if not np.isnan(ratio[i]):
                self.assertTrue(ratio[i] < 1.0)

    def test_gaussian_background(self): 
        gen = GaussianBackgroundGenerator(self.bkgd)
        dev = next(gen)
        self.assertEqual(dev.size, self.bkgd.size)
        self.assertCountEqual(dev.lo_edges, self.bkgd.lo_edges)
        self.assertCountEqual(dev.hi_edges, self.bkgd.hi_edges)
        self.assertCountEqual(dev.exposure, self.bkgd.exposure)
        devs = [next(gen) for i in range(10)]

    def test_variable_poisson_background(self): 
        gen = VariableGaussianBackground(self.bkgd)
        dev = next(gen)
        self.assertEqual(dev.size, self.bkgd.size)
        self.assertCountEqual(dev.lo_edges, self.bkgd.lo_edges)
        self.assertCountEqual(dev.hi_edges, self.bkgd.hi_edges)
        self.assertCountEqual(dev.exposure, self.bkgd.exposure)
        devs = [next(gen) for i in range(10)]
        
        gen.amp = 10.0
        dev = next(gen)
        ratio = dev.counts/devs[0].counts
        for i in range(self.bkgd.size):
            if not np.isnan(ratio[i]):
                self.assertTrue(ratio[i] > 1.0)

        gen.amp = 0.1
        dev = next(gen)
        ratio = dev.counts/devs[0].counts
        for i in range(self.bkgd.size):
            if not np.isnan(ratio[i]):
                self.assertTrue(ratio[i] < 1.0)
    
    def test_source_spectrum(self):
        gen = SourceSpectrumGenerator(self.rsp, self.fxn, self.params, self.exposure)
        dev = next(gen)
        self.assertEqual(dev.size, self.rsp.numchans)
        self.assertCountEqual(dev.lo_edges, self.rsp.ebounds['E_MIN'])
        self.assertCountEqual(dev.hi_edges, self.rsp.ebounds['E_MAX'])
        self.assertEqual(dev.exposure[0], self.exposure)
        devs = [next(gen) for i in range(10)]

    def test_variable_source_spectrum(self):
        gen = VariableSourceSpectrumGenerator(self.rsp, self.fxn, self.params, 
                                              self.exposure)
        dev = next(gen)
        self.assertEqual(dev.size, self.rsp.numchans)
        self.assertCountEqual(dev.lo_edges, self.rsp.ebounds['E_MIN'])
        self.assertCountEqual(dev.hi_edges, self.rsp.ebounds['E_MAX'])
        self.assertEqual(dev.exposure[0], self.exposure)
        devs = [next(gen) for i in range(10)]
        
        gen.amp = 10.0
        dev = next(gen)
        ratio = dev.counts/devs[0].counts
        for i in range(self.bkgd.size):
            if not np.isnan(ratio[i]):
                self.assertTrue(ratio[i] > 1.0)
    
    def test_event_spectrum(self):
        spec_gen = PoissonBackgroundGenerator(self.bkgd)
        event_gen = EventSpectrumGenerator(next(spec_gen).counts, 0.001)
        dev_times, dev_chans = next(event_gen)
        self.assertEqual(dev_times.size, dev_chans.size)
        
        event_gen.spectrum = next(spec_gen).counts
        dev_times, dev_chans = next(event_gen)
        self.assertEqual(dev_times.size, dev_chans.size)


class TestPhaSimulator(TestCase):
    bkgd = BAK.open(os.path.join(data_dir, 'glg_tte_n0_bn160509374_xspec_v00.bak'))
    bkgd = bkgd.data
    rsp = RSP.open(os.path.join(data_dir, 'glg_cspec_n4_bn120415958_v00.rsp2'))
    fxn = Band()
    params = (1.0, 300.0, -1.0, -2.5)
    exposure = 2.048
    
    def test_run(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        
        sim = pha_sims.simulate_background(1)
        self.assertEqual(sim[0].size, self.bkgd.size)
        sim = pha_sims.simulate_source(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
        sim = pha_sims.simulate_sum(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
            
        sims = pha_sims.simulate_sum(10)
        self.assertEqual(len(sims), 10)
    
    def test_set_rsp(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        pha_sims.set_rsp(self.rsp.extract_drm(5)) 
        sim = pha_sims.simulate_background(1)
        self.assertEqual(sim[0].size, self.bkgd.size)
        sim = pha_sims.simulate_source(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
        sim = pha_sims.simulate_sum(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)

    def test_set_background(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        pha_sims.set_background(self.bkgd, 'Poisson') 
        sim = pha_sims.simulate_background(1)
        self.assertEqual(sim[0].size, self.bkgd.size)
        sim = pha_sims.simulate_source(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
        sim = pha_sims.simulate_sum(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)

    def test_set_source(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        pha_sims.set_source(self.fxn, self.params, 10.0) 
        sim = pha_sims.simulate_background(1)
        self.assertEqual(sim[0].size, self.bkgd.size)
        sim = pha_sims.simulate_source(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
        sim = pha_sims.simulate_sum(1)
        self.assertEqual(sim[0].size, self.rsp.numchans)
    
    def test_to_bak(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        baks = pha_sims.to_bak(5)
        [self.assertIsInstance(bak, BAK) for bak in baks]
        baks = pha_sims.to_bak(5, tstart=5.0, tstop=8.0)
        [self.assertIsInstance(bak, BAK) for bak in baks]

    def test_to_pha(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        phas = pha_sims.to_pha(5)
        [self.assertIsInstance(pha, PHA) for pha in phas]
        phas = pha_sims.to_pha(5, tstart=5.0, tstop=8.0)
        [self.assertIsInstance(pha, PHA) for pha in phas]

    def test_to_phaii(self):
        pha_sims = PhaSimulator(self.rsp, self.fxn, self.params, self.exposure,
                                self.bkgd, 'Gaussian')
        phaii = pha_sims.to_phaii(5)
        self.assertIsInstance(phaii, PHAII)
        phaii = pha_sims.to_phaii(5, bin_width=3.0)
        self.assertIsInstance(phaii, PHAII)
     

class TestTteSourceSimulator(TestCase):
    rsp = RSP.open(os.path.join(data_dir, 'glg_cspec_n4_bn120415958_v00.rsp2'))
    spec_fxn = Band()
    spec_params = (0.1, 300.0, -1.0, -2.5)
    time_fxn = tophat
    time_params = (0.05, 0.0, 1.0)
    
    def test_run(self):
        src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params)
        tte = src_sims.to_tte(-5.0, 5.0, trigtime=0.0)
        self.assertIsInstance(tte, TTE)
        tte = src_sims.to_tte(-1.0, 1.0)
        self.assertIsInstance(tte, TTE)
    
    def test_set_response(self):
        src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params)
        src_sims.set_response(self.rsp.extract_drm(5))   
        tte = src_sims.to_tte(-1.0, 1.0)
        self.assertIsInstance(tte, TTE)

    def test_set_spectrum(self):
        src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params)
        src_sims.set_spectrum(self.spec_fxn, (0.1, 100., -0.5, -3.0))   
        tte = src_sims.to_tte(-1.0, 1.0)
        self.assertIsInstance(tte, TTE)

    def test_time_profile(self):
        src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params)
        src_sims.set_time_profile(tophat, (0.1, 0.0, 2.0))   
        tte = src_sims.to_tte(-1.0, 1.0)
        self.assertIsInstance(tte, TTE)
    
    def test_errors(self):
        with  self.assertRaises(ValueError):
            src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params, sample_period=0.0)
        
        with  self.assertRaises(ValueError):
            src_sims = TteSourceSimulator(self.rsp, self.spec_fxn, self.spec_params,
                                      tophat, self.time_params, deadtime=-1e-4)

        

class TestTteBackgroundSimulator(TestCase):
    bkgd = BAK.open(os.path.join(data_dir, 'glg_tte_n0_bn160509374_xspec_v00.bak'))
    bkgd = bkgd.data
    time_fxn = linear
    time_params = (1.0, -0.1)
    
    def test_run(self):
        src_sims = TteBackgroundSimulator(self.bkgd, 'Gaussian', linear, 
                                          self.time_params)
        tte = src_sims.to_tte(-5.0, 5.0)
        self.assertIsInstance(tte, TTE)
        tte = src_sims.to_tte(-1.0, 1.0, trigtime=0)
        self.assertIsInstance(tte, TTE)
    
    def test_set_background(self):
        src_sims = TteBackgroundSimulator(self.bkgd, 'Gaussian', linear, 
                                          self.time_params)
        src_sims.set_background(self.bkgd, 'Poisson')   
        tte = src_sims.to_tte(-1.0, 1.0)
        self.assertIsInstance(tte, TTE)

    def test_errors(self):
        with  self.assertRaises(ValueError):
            src_sims = TteBackgroundSimulator(self.bkgd, 'Gaussian', linear, 
                                              self.time_params, sample_period=0.0)
        
        with  self.assertRaises(ValueError):
            src_sims = TteBackgroundSimulator(self.bkgd, 'Gaussian', linear, 
                                              self.time_params, deadtime=-1e-4)


class TestTimeProfiles(TestCase):
    times = np.array([-10.0, 0.0, 10.0])
    
    def test_tophat(self):
        params = (1.0, 0.0, 20.0)
        y = tophat(self.times, *params)
        self.assertCountEqual(y, np.array([0.0, 1.0, 1.0]))
    
    def test_norris(self):
        params = (1.0, -1.0, 0.1, 2.0)
        y = norris(self.times, *params)
        true = np.array((0.0, 0.858, 0.006))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]
    
    def test_constant(self):
        params = (1.0,)
        y = constant(self.times, *params)
        self.assertCountEqual(y, np.array([1.0, 1.0, 1.0]))

    def test_linear(self):
        params = (1.0, -2.0,)
        y = linear(self.times, *params)
        self.assertCountEqual(y, np.array([21.0, 1.0, -19.0]))

    def test_quadratic(self):
        params = (1.0, -2.0, 2.0)
        y = quadratic(self.times, *params)
        self.assertCountEqual(y, np.array([221.0, 1.0, 181.0]))
        

if __name__ == '__main__':
    unittest.main()
      
        
        