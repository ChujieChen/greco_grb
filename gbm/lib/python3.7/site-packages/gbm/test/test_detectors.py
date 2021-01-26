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
from gbm.detectors import *
from copy import copy


class TestDetectors(unittest.TestCase):

    expected_nai = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb']
    expected_bgo = ['b0', 'b1']
    expected_all = expected_nai + expected_bgo

    def test_detector_list(self):
        detectors = copy(self.expected_all)
        for d in Detector:
            if d.short_name in detectors:
                detectors.remove(d.short_name)
            else:
                self.fail("Detector was not expected")
        if detectors:
            self.fail("All of the detectors weren't removed")

    def test_is_nai(self):
        for d in self.expected_nai:
            self.assertTrue(Detector.is_nai(Detector.from_str(d)))
        for d in self.expected_bgo:
            self.assertFalse(Detector.is_nai(Detector.from_str(d)))

    def test_is_bgo(self):
        for d in self.expected_bgo:
            self.assertTrue(Detector.is_bgo(Detector.from_str(d)))
        for d in self.expected_nai:
            self.assertFalse(Detector.is_bgo(Detector.from_str(d)))

    def test_nai(self):
        l = Detector.nai()
        self.assertEqual(len(l), len(self.expected_nai))
        for d in self.expected_nai:
            self.assertIsNotNone(Detector.from_str(d))

    def test_bgo(self):
        l = Detector.bgo()
        self.assertEqual(len(l), len(self.expected_bgo))
        for d in self.expected_bgo:
            self.assertIsNotNone(Detector.from_str(d))

    def test_from_number(self):
        num = 0
        for d in self.expected_all:
            self.assertEqual(d, Detector.from_num(num).short_name)
            num += 1

    def test_meta_data(self):
        expected_data = [
            ('N0', 'NAI_00', 0, 45.89, 20.58),
            ('N1', 'NAI_01', 1, 45.11, 45.31),
            ('N2', 'NAI_02', 2, 58.44, 90.21),
            ('N3', 'NAI_03', 3, 314.87, 45.24),
            ('N4', 'NAI_04', 4, 303.15, 90.27),
            ('N5', 'NAI_05', 5, 3.35, 89.79),
            ('N6', 'NAI_06', 6, 224.93, 20.43),
            ('N7', 'NAI_07', 7, 224.62, 46.18),
            ('N8', 'NAI_08', 8, 236.61, 89.97),
            ('N9', 'NAI_09', 9, 135.19, 45.55),
            ('NA', 'NAI_10', 10, 123.73, 90.42),
            ('NB', 'NAI_11', 11, 183.74, 90.32),
            ('B0', 'BGO_00', 12, 0.00, 90.00),
            ('B1', 'BGO_01', 13, 180.00, 90.00),
        ]
        for d in expected_data:
            det = Detector.from_num(d[2])
            self.assertEqual(det.name,  d[0])
            self.assertEqual(det.long_name, d[1])
            self.assertEqual(det.azimuth, d[3])
            self.assertEqual(det.zenith, d[4])
            self.assertEqual(det.pointing, (d[3], d[4]))
            self.assertEqual(det.__repr__(), "Detector(\"{}\", \"{}\", {})".format(d[0], d[1], d[2]))

    def test_invalid_str(self):
        self.assertIsNone(Detector.from_str("FAKE"))

    def test_invalid_num(self):
        self.assertIsNone(Detector.from_num(20))

if __name__ == '__main__':
    unittest.main()