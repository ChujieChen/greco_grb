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

from gbm.data.poshist import PosHist

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestPosHist(TestCase):
    filename = os.path.join(data_dir, 'glg_poshist_all_170101_v00.fit')
    tstart = 504921540.740104
    tstop = 505008061.340078
    test_time = 504931642.739272 #trigger 170101.116
    def test_attributes(self):
        p = PosHist.open(self.filename)
        self.assertCountEqual(p.headers.keys(), ['PRIMARY', 'GLAST POS HIST'])
        self.assertCountEqual(p.time_range, (self.tstart, self.tstop))
        self.assertEqual(len(p.gti), 10)
        #mark TODO: Need to figure out what to do for daily poshist flags
        # The flags are incorrect in the file: The Sun is never visible
        
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
        
        p = PosHist.open(self.filename)
        
        eic = p.get_eic(self.test_time)/1000.0
        [self.assertAlmostEqual(eic[i], test_eic[i], delta=1.0) for i in range(3)]
        quat = p.get_quaternions(self.test_time)
        [self.assertAlmostEqual(quat[i], test_quat[i], places=4) for i in range(4)]
        
        lat = p.get_latitude(self.test_time)
        lon = p.get_longitude(self.test_time)
        alt = p.get_altitude(self.test_time)/1000.0
        self.assertAlmostEqual(lat, test_lat, delta=0.1)
        self.assertAlmostEqual(lon, test_lon, delta=0.1)
        self.assertAlmostEqual(alt, test_alt, delta=5.0)
        
        vel = p.get_velocity(self.test_time)
        [self.assertAlmostEqual(vel[i]/1000.0, test_vel[i], delta=0.01) for i in range(3)]
        angvel = p.get_angular_velocity(self.test_time)
        [self.assertAlmostEqual(angvel[i], test_angvel[i], places=5) for i in range(3)]
        
        geo_radec = p.get_geocenter_radec(self.test_time)
        self.assertAlmostEqual(geo_radec[0], test_geoloc[0], delta=0.1)
        self.assertAlmostEqual(geo_radec[1], test_geoloc[1], delta=0.1)
        geo_radius = p.get_earth_radius(self.test_time)
        self.assertAlmostEqual(geo_radius, test_earth_radius, delta=0.01)
        
        sun_visible = p.get_sun_visibility(self.test_time)
        #self.assertEqual(sun_visible, False)
        in_saa = p.get_saa_passage(self.test_time)
        self.assertEqual(in_saa, False)
        ml = p.get_mcilwain_l(self.test_time)
        self.assertAlmostEqual(ml, test_mcilwain, delta=0.01)
    
    def test_coordinate_conversions(self):
        test_ra, test_dec = 70.64, -1.58
        test_az, test_zen = 138.0, 65.0
        test_n0_ra, test_n0_dec = 30.908, 57.679
        test_n0_angle = 67.42
        
        p = PosHist.open(self.filename)
        
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
            p = PosHist.open('roland.fit')
        
        p = PosHist.open(self.filename)
        with self.assertRaises(ValueError):
            p.get_eic(1.0)
            p.detector_pointing('n11', 1.0)
            p.detector_angle(0.0, 0.0, 'a1', 1.0)

if __name__ == '__main__':
    unittest.main()