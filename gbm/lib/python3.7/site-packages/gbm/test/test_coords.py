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
from gbm.coords import *


class TestCoords(unittest.TestCase):

    def test_haversine(self):
        angle = haversine(10.0, 0.0, 20.0, 0.0)
        self.assertEqual(angle, 10.0)
        angle = haversine(0.0, -10.0, 0.0, -20.0)
        self.assertEqual(angle, 10.0)
        angle = np.rad2deg(haversine(np.deg2rad(10.0), np.deg2rad(0.0), \
                            np.deg2rad(20.0), np.deg2rad(0.0), deg=False))
        self.assertEqual(angle, 10.0)
        
    def test_sphere2cartesian(self):
        true_cart = np.array([0.25, 0.4330127019, 0.8660254038])
        cart = azzen_to_cartesian(60.0, 30.0)
        for i in range(3):
            self.assertAlmostEqual(cart[i], true_cart[i])
        
        cart = radec_to_cartesian(60.0, 90.0-30.0)
        for i in range(3):
            self.assertAlmostEqual(cart[i], true_cart[i])

    def test_quaternions(self):
        quat = np.array([1.0, 1.0, 1.0, 1.0])
        true_conj = np.array([-1.0, -1.0, -1.0, 1.0])
        conj = quaternion_conj(quat)
        for i in range(4):
            self.assertEqual(conj[i], true_conj[i])
        
        true_inverse = np.array([-0.25, -0.25, -0.25, 0.25])
        inverse = quaternion_inv(quat)
        for i in range(4):
            self.assertEqual(inverse[i], true_inverse[i])
        
        quat2 = np.array([-1.0, -1.0, -1.0, -1.0])
        true_product = np.array([-2.0, -2.0, -2.0, 2.0])
        product = quaternion_prod(quat, quat2)
        for i in range(4):
            self.assertEqual(product[i], true_product[i])
        
        true_dcm = np.array([[0,0,1], [1,0,0], [0,1,0]])
        dcm = spacecraft_direction_cosines(quat)
        for i in range(3):
            for j in range(3):
                self.assertEqual(dcm[i,j], true_dcm[i,j])
    
    def test_geocenter(self):
        coord = np.array([-3227.0, 6092.0, 471.0])
        true_geocenter = np.array([297.911, -3.908])
        geocenter = geocenter_in_radec(coord)
        for i in range(2):
            self.assertAlmostEqual(geocenter[i], true_geocenter[i], 3)
    
    def test_fermi2radec(self):
        quat = np.array([1.0, 1.0, 1.0, 1.0])
        az = 180.0
        zen = 0.0
        true_loc = np.array([0.0, 0.0])
        loc = spacecraft_to_radec(az, zen, quat)
        for i in range(2):
            self.assertEqual(loc[i], true_loc[i])
                
    def test_radec2fermi(self):
        quat = np.array([1.0, 1.0, 1.0, 1.0])
        ra = 180.0
        dec = 0.0
        true_loc = np.array([0.0, 180.0])
        loc = radec_to_spacecraft(ra, dec, quat)
        for i in range(2):
            self.assertEqual(loc[i], true_loc[i])
    
    def test_latitude(self):
        coord = np.array([-3227.0, 6092.0, 471.0])*1000.0
        true_lat = 3.91
        true_alt = 538.972*1000.0
        lat, alt = latitude_from_geocentric_coords_simple(coord)
        self.assertAlmostEqual(lat, true_lat, 2)
        self.assertAlmostEqual(np.round(alt), true_alt, 2)
    
        true_lat = 3.93
        true_alt = 531.944*1000.0
        lat, alt = latitude_from_geocentric_coords_complex(coord)
        self.assertAlmostEqual(lat, true_lat, 2)
        self.assertAlmostEqual(np.round(alt), true_alt, 2)
    
    def test_longitude(self):
        coord = np.array([-3227.0, 6092.0, 471.0])
        met = 524666471.0
        true_lon = 321.552
        lon = longitude_from_geocentric_coords(coord, met)
        self.assertAlmostEqual(lon, true_lon, 3)
        
        lon = longitude_from_geocentric_coords(coord, met, ut1=True)
        self.assertAlmostEqual(lon, true_lon, 2)
        
    def test_mcilwainl(self):
        lon = 321.552
        lat = 3.93
        true_mcilwainl = 1.12
        mcilwainl = calc_mcilwain_l(lat, lon)
        self.assertAlmostEqual(mcilwainl, true_mcilwainl, 2)
        
    def test_sun_loc(self):
        met = 524666471.0
        true_ra = 146.85
        true_dec = 13.34
        ra, dec = get_sun_loc(met)
        self.assertAlmostEqual(ra, true_ra, 2)
        self.assertAlmostEqual(dec, true_dec, 2)

if __name__ == '__main__':
    unittest.main()
      
