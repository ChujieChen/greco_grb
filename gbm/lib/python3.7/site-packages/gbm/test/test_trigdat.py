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

from gbm.data.trigdat import Trigdat
from gbm.data.tcat import Tcat
#from gbm.binning.unbinned import bin_by_time
#from gbm.binning.binned import combine_by_factor

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestTrigdat(TestCase):
    filename = os.path.join(data_dir, 'glg_trigdat_all_bn170101116_v01.fit')
    trigtime = 504931642.867272
    header_names = ['PRIMARY', 'TRIGRATE', 'BCKRATES', 'OB_CALC', 'MAXRATES', 'EVNTRATE']
    num_maxrates = 3
    time_range = (504931505.008880, 504932119.419874)
    test_time = 504931642.739272 #trigger 170101.116
    
    def test_attributes(self):
        t = Trigdat.open(self.filename)
        self.assertEqual(t.trigtime, self.trigtime)
        self.assertEqual(t.is_gbm_file, True)
        self.assertEqual(t.id, '170101116')
        self.assertEqual(t.filename, os.path.basename(self.filename))
        self.assertEqual(t.is_trigger, True)
        self.assertEqual(t.detector, 'all')
        self.assertEqual(t.datatype, 'TRIGDAT')
        self.assertCountEqual(t.headers.keys(), self.header_names)
        self.assertEqual(t.num_maxrates, self.num_maxrates)
        self.assertAlmostEqual(t.time_range[0], self.time_range[0], places=6)
        self.assertAlmostEqual(t.time_range[1], self.time_range[1], places=6)
        self.assertEqual(t.maxrates[0].numchans, 8)
        self.assertEqual(t.maxrates[0].numdets, 14)
        self.assertEqual(t.maxrates[0].time_range, (504931642.227346, 504931646.323346))
        self.assertAlmostEqual(t.maxrates[0].timescale, 4.096)
        self.assertEqual(t.backrates.numchans, 8)
        self.assertEqual(t.backrates.numdets, 14)
        self.assertEqual(t.backrates.time_range, (504931609.913434, 504931642.68143404))
        loc = t.fsw_locations[0]
        self.assertEqual(loc.time, 504931644.275308)
        self.assertAlmostEqual(loc.location[0], 64.25, places=4)
        self.assertAlmostEqual(loc.location[1], 8.55, places=4)
        self.assertAlmostEqual(loc.location[2], 8.016666, places=4)
        self.assertEqual(loc.top_classification[0], 'GRB')
        self.assertAlmostEqual(loc.top_classification[1], 0.854902, places=4)
        self.assertEqual(loc.next_classification[0], 'GROJ422')
        self.assertAlmostEqual(loc.next_classification[1], 0.1019608, places=4)
        self.assertAlmostEqual(loc.intensity, 277.1126, places=4)
        self.assertAlmostEqual(loc.hardness_ratio, 0.5011158, places=4)
        self.assertAlmostEqual(loc.fluence, 0, places=4)
        self.assertAlmostEqual(loc.significance, 19.9, places=4)
        self.assertAlmostEqual(loc.timescale, 2.048, places=4)
        self.assertEqual(loc.spectrum, 'soft')
        self.assertEqual(loc.location_sc, (130, 55))
        self.assertEqual(t.triggered_detectors, ['n9', 'na', 'nb'])
    
    def test_get_maxrates(self):
        t = Trigdat.open(self.filename)
        self.assertEqual(t.get_maxrates(1), t.maxrates[1])        
        self.assertEqual(t.get_fsw_locations(1), t.fsw_locations[1])
        
    def test_to_ctime(self):
        t = Trigdat.open(self.filename)
        c = t.to_ctime('n5')
        self.assertEqual(c.detector, 'n5')
        self.assertEqual(c.datatype, 'CTIME')
        self.assertAlmostEqual(c.time_range[0], t.time_range[0]-t.trigtime)
        self.assertAlmostEqual(c.time_range[1], t.time_range[1]-t.trigtime)
        self.assertEqual(c.trigtime, t.trigtime)
        self.assertEqual(c.numchans, 8)
    
    def test_sum_detectors(self):
        t = Trigdat.open(self.filename)
        c1 = t.to_ctime('n5')
        c2 = t.sum_detectors(['n5', 'n5'])
        self.assertCountEqual(c1.data.counts.flatten()*2, 
                              c2.data.counts.flatten())
    
    def test_interpolators(self):
        test_eic = np.array([-1572., 6370., 2164.])
        test_quat = np.array([0.2229096, 0.06231983, 0.5392869, -0.8096896])
        test_lat, test_lon = 18.23, 321.02
        test_alt = np.sqrt(np.sum(test_eic**2))-6371.0 # simple altitude calc
        test_vel = np.array([-6820.45, -2465.3848, 2270.4202])/1000.0
        test_angvel = np.array([0.000706, 6.29781e-06, -0.000714])
        test_geoloc = np.array([283.85, -18.25])
        test_earth_radius = 67.335
        test_mcilwain = 1.22
        
        p = Trigdat.open(self.filename)
        
        eic = p.get_eic(self.test_time)/1000.0
        [self.assertAlmostEqual(eic[i], test_eic[i], delta=5.0) for i in range(3)]
        quat = p.get_quaternions(self.test_time)
        [self.assertAlmostEqual(quat[i], test_quat[i], places=2) for i in range(4)]
        
        lat = p.get_latitude(self.test_time)
        lon = p.get_longitude(self.test_time)
        alt = p.get_altitude(self.test_time)/1000.0
        self.assertAlmostEqual(lat, test_lat, delta=1.0)
        self.assertAlmostEqual(lon, test_lon, delta=1.0)
        self.assertAlmostEqual(alt, test_alt, delta=5.0)
        
        vel = p.get_velocity(self.test_time)
        [self.assertAlmostEqual(vel[i]/1000.0, test_vel[i], delta=1.0) for i in range(3)]
        angvel = p.get_angular_velocity(self.test_time)
        [self.assertAlmostEqual(angvel[i], test_angvel[i], places=3) for i in range(3)]
        
        geo_radec = p.get_geocenter_radec(self.test_time)
        self.assertAlmostEqual(geo_radec[0], test_geoloc[0], delta=0.1)
        self.assertAlmostEqual(geo_radec[1], test_geoloc[1], delta=0.1)
        geo_radius = p.get_earth_radius(self.test_time)
        self.assertAlmostEqual(geo_radius, test_earth_radius, delta=0.01)
        
        sun_visible = p.get_sun_visibility(self.test_time)
        self.assertEqual(sun_visible, False)
        in_saa = p.get_saa_passage(self.test_time)
        self.assertEqual(in_saa, False)
        ml = p.get_mcilwain_l(self.test_time)
        self.assertAlmostEqual(ml, test_mcilwain, delta=0.01)

    def test_coordinate_conversions(self):
        test_ra, test_dec = 70.64, -1.58
        test_az, test_zen = 138.0, 65.0
        test_n0_ra, test_n0_dec = 30.908, 57.679
        test_n0_angle = 67.42
        
        p = Trigdat.open(self.filename)
        
        az, zen = p.to_fermi_frame(test_ra, test_dec, self.test_time)
        self.assertAlmostEqual(az, test_az, delta=1.0)
        self.assertAlmostEqual(zen, test_zen, delta=1.0)
        
        ra, dec = p.to_equatorial(test_az, test_zen, self.test_time)
        self.assertAlmostEqual(ra, test_ra, delta=1.0)
        self.assertAlmostEqual(dec, test_dec, delta=1.0)
        
        loc_vis = p.location_visible(test_ra, test_dec, self.test_time)
        self.assertEqual(loc_vis, True)
        
        ra, dec = p.detector_pointing('n0', self.test_time)
        self.assertAlmostEqual(ra, test_n0_ra, delta=0.5)
        self.assertAlmostEqual(dec, test_n0_dec, delta=0.5)
        angle = p.detector_angle(test_ra, test_dec, 'n0', self.test_time)
        self.assertAlmostEqual(angle, test_n0_angle, delta=0.5)

    def test_errors(self):
        with  self.assertRaises(IOError):
            t = Trigdat.open('deschain.fit')
        
        t = Trigdat.open(self.filename)
        with self.assertRaises(ValueError):
            t.to_ctime('b4')
            t.to_ctime('n0', timescale=123)
