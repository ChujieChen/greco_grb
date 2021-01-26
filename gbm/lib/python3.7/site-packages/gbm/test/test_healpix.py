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

from gbm.data.localization import *

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


class TestHealPix(TestCase):
    filename = os.path.join(data_dir, 'glg_healpix_all_bn190915240_v00.fit')
    trigtime = 590219102.911008
    npix = 196608
    nside = 128
    pixel_area = 0.2098
    sun = (172.5011935415178, 3.23797213866954)
    geo = (319.8312390218318, 17.40612934717674)
    geo_radius = 67.2950460311874
    n0 = (146.5959532829778, 36.96759511828569)
    scpos = [-5039500.,  4254000., -2067500.]
    quat = [-0.223915,  0.447149,  0.860062, -0.101055]
    centroid = (48.8671875, 4.181528273111476)
    
    def test_attributes(self):
        h = GbmHealPix.open(self.filename)
        self.assertEqual(h.trigtime, self.trigtime)
        self.assertEqual(h.npix, self.npix)
        self.assertEqual(h.nside, self.nside)
        self.assertAlmostEqual(h.pixel_area, self.pixel_area, places=4)
        self.assertEqual(h.sun_location, self.sun)
        self.assertEqual(h.geo_location, self.geo)
        self.assertEqual(h.geo_radius, self.geo_radius)
        self.assertEqual(h.n0_pointing, self.n0)
        self.assertEqual(h.scpos.tolist(), self.scpos)
        self.assertEqual(h.quaternion.tolist(), self.quat)
    
    def test_geo(self):
        h = GbmHealPix.open(self.filename)
        self.assertAlmostEqual(h.geo_probability, 0.0, places=4)
        self.assertAlmostEqual(h.observable_fraction(h), 1.0, places=4)
        h2 = h.remove_earth(h)
        self.assertEqual(h2.geo_probability, 0.0)
        self.assertEqual(h2._prob.sum(), 1.0)       
        self.assertAlmostEqual(h2.observable_fraction(h2), 1.0, places=4)
    
    def test_assoc_prob(self):
        h = GbmHealPix.open(self.filename)
        self.assertAlmostEqual(h.source_probability(*h.centroid), 0.997, places=3)
        self.assertAlmostEqual(h.region_probability(h), 0.995, places=3)
        
        self.assertEqual(h.source_probability(*h.centroid, prior=0.0), 0.0)
        self.assertEqual(h.source_probability(*h.centroid, prior=1.0), 1.0)
        self.assertEqual(h.region_probability(h, prior=0.0), 0.0)
        self.assertEqual(h.region_probability(h, prior=1.0), 1.0)
        
        g = GbmHealPix.from_gaussian(100.0, 60.0, 5.0, nside=128)
        self.assertAlmostEqual(h.region_probability(g), 0.0, places=6)
        self.assertEqual(h.source_probability(100.0, 60.0), 0.0)
        
        
    def test_annulus(self):
        h = GbmHealPix.from_annulus(300.0, 10.0, 70.0, 10.0, nside=128)
        self.assertEqual(h.nside, self.nside)
        self.assertEqual(h.npix, self.npix)
        
    def test_multiply(self):
        h1 = GbmHealPix.open(self.filename)
        h2 = GbmHealPix.from_annulus(300.0, 10.0, 70.0, 10.0)
        h = GbmHealPix.multiply(h1, h2)
        self.assertEqual(h.nside, self.nside)
        self.assertEqual(h.npix, self.npix)

class TestChi2Grid(TestCase):
    filename = os.path.join(data_dir, 'chi2grid_bn190531568_v00.dat')
    trigtime = 581002688.0
    scpos = [5761500., -3302750., 1907250.]
    quat = [-0.592056, -0.192717, -0.517305, -0.587134]
    numpts = 41168
    
    def test_attributes(self):
        c = Chi2Grid.open(self.filename)
        self.assertEqual(c.numpts, self.numpts)
        self.assertEqual(c.azimuth.size, self.numpts)
        self.assertEqual(c.zenith.size, self.numpts)
        self.assertEqual(c.ra.size, self.numpts)
        self.assertEqual(c.dec.size, self.numpts)
        self.assertEqual(c.chisq.size, self.numpts)
        self.assertEqual(c.significance.size, self.numpts)
        
    def test_from_data(self):
        c = Chi2Grid.open(self.filename)
        c2 = Chi2Grid.from_data(c.azimuth, c.zenith, c.ra, c.dec, c.chisq)
        self.assertCountEqual(c.azimuth, c2.azimuth)
        self.assertCountEqual(c.zenith, c2.zenith)
        self.assertCountEqual(c.ra, c2.ra)
        self.assertCountEqual(c.dec, c2.dec)
        self.assertCountEqual(c.chisq, c2.chisq)
        self.assertCountEqual(c.significance, c2.significance)
    
    def test_healpix_from_chi2grid(self):
        c = Chi2Grid.open(self.filename)
        c.quaternion = self.quat
        c.scpos = self.scpos
        c.trigtime = self.trigtime
        
        h = GbmHealPix.from_chi2grid(c, nside=512)
        self.assertEqual(h.quaternion.tolist(), self.quat)
        self.assertEqual(h.scpos.tolist(), self.scpos)
        self.assertEqual(h.trigtime, self.trigtime)
        self.assertAlmostEqual(h.centroid[0], c.ra[c.chisq.argmin()], places=-1)
        self.assertAlmostEqual(h.centroid[1], c.dec[c.chisq.argmin()], places=-1)

class TestSystematics(TestCase):
    
    def test_gbuts(self):
        sig, frac = GBUTS_Model_O3()
        self.assertEqual(frac, [1.0])
        self.assertAlmostEqual(sig[0], 0.047, places=3)

    def test_hitl(self):
        sig, frac = HitL_Model(50.0)
        self.assertEqual(frac, [0.918])
        self.assertAlmostEqual(sig[0], 0.073, places=3)
        self.assertAlmostEqual(sig[1], 0.267, places=3)

        sig, frac = HitL_Model(100.0)
        self.assertEqual(frac, [0.884])
        self.assertAlmostEqual(sig[0], 0.040, places=3)
        self.assertAlmostEqual(sig[1], 0.230, places=3)

    def test_ground(self):
        sig, frac = GA_Model()
        self.assertEqual(frac, [0.804])
        self.assertAlmostEqual(sig[0], 0.065, places=3)
        self.assertAlmostEqual(sig[1], 0.239, places=3)

    def test_roboba(self):
        sig, frac = RoboBA_Function('long')
        self.assertEqual(frac, [0.579])
        self.assertAlmostEqual(sig[0], 0.032, places=3)
        self.assertAlmostEqual(sig[1], 0.072, places=3)

        sig, frac = RoboBA_Function('short')
        self.assertEqual(frac, [0.390])
        self.assertAlmostEqual(sig[0], 0.045, places=3)
        self.assertAlmostEqual(sig[1], 0.077, places=3)

    def test_untargeted(self):
        sig, frac = Untargeted_Search_Model()
        self.assertEqual(frac, [1.0])
        self.assertAlmostEqual(sig[0], 0.097, places=3)
                

if __name__ == '__main__':
    unittest.main()
      
