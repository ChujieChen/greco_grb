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

from gbm.data.phaii import Cspec
from gbm.data.primitives import TimeBins
from gbm.data.collection import DataCollection, GbmDetectorCollection

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 
                        '160509374')

class TestDataCollection(TestCase):
    b0 = Cspec.open(os.path.join(data_dir, 'glg_cspec_b0_bn160509374_v01.pha'))
    n0 = Cspec.open(os.path.join(data_dir, 'glg_cspec_n0_bn160509374_v01.pha'))
    n1 = Cspec.open(os.path.join(data_dir, 'glg_cspec_n1_bn160509374_v01.pha'))
    collection = DataCollection.from_list([b0, n0, n1])
    

    def test_attributes(self):
        self.assertEqual(len(self.collection), 3)
        self.assertCountEqual(self.collection.items, [self.b0.filename, 
                                                      self.n0.filename,
                                                      self.n1.filename])
        self.assertEqual(self.collection.types, Cspec)
    
    def test_get_item(self):
        item = self.collection.get_item(self.collection.items[0])
        self.assertEqual(item, self.b0)
    
    def test_remove_and_include(self):
        # remove
        self.collection.remove(self.collection.items[2])
        self.assertEqual(len(self.collection), 2)
        self.assertCountEqual(self.collection.items, [self.b0.filename, 
                                                      self.n0.filename])
        # include
        self.collection.include(self.n1)
        self.test_attributes()
    
    def test_to_list(self):
        thelist = self.collection.to_list()
        self.assertCountEqual(thelist, [self.b0, self.n0, self.n1])

    def test_item_attributes(self):
        numchans = self.collection.numchans()
        self.assertCountEqual(numchans, [128, 128, 128])
        erange = self.collection.energy_range()
        elo = (113.00731, 4.5702357, 4.089358)
        ehi = (50000.0, 2000.0, 2000.0)
        [self.assertAlmostEqual(erange[i][0], elo[i], places=3) for i in range(3)]
        [self.assertAlmostEqual(erange[i][1], ehi[i], places=3) for i in range(3)]

    def test_item_methods(self):
        exposure = self.collection.get_exposure()
        test_exp = [7967., 7979., 7978.]
        [self.assertAlmostEqual(exposure[i], test_exp[i], places=0) for i in range(3)]

        lcs = self.collection.to_lightcurve()
        [self.assertIsInstance(lc, TimeBins) for lc in lcs]


class TestGbmDetectorCollection(TestCase):
    b0 = Cspec.open(os.path.join(data_dir, 'glg_cspec_b0_bn160509374_v01.pha'))
    n0 = Cspec.open(os.path.join(data_dir, 'glg_cspec_n0_bn160509374_v01.pha'))
    n1 = Cspec.open(os.path.join(data_dir, 'glg_cspec_n1_bn160509374_v01.pha'))
    collection = GbmDetectorCollection.from_list([b0, n0, n1])
    
    def test_attributes(self):
        self.assertEqual(len(self.collection), 3)
        self.assertCountEqual(self.collection.items, [self.b0.filename, 
                                                      self.n0.filename,
                                                      self.n1.filename])
        self.assertEqual(self.collection.types, Cspec)

    def test_remove_and_include(self):
        # remove
        self.collection.remove(self.collection.items[2])
        self.assertEqual(len(self.collection), 2)
        self.assertCountEqual(self.collection.items, [self.b0.filename, 
                                                      self.n0.filename])
        # include
        self.collection.include(self.n1, 'n1')
        self.test_attributes()
    
    def test_item_methods(self):
        slices = self.collection.slice_time(nai_args=((0.0, 10.0),), 
                                            bgo_args=((0.0, 100.0),))
        test_exp = [102., 12., 12.]
        [self.assertAlmostEqual(slices[i].get_exposure(), test_exp[i], places=0) \
         for i in range(3)]
        
        specs = self.collection.to_spectrum(nai_kwargs={'energy_range':(8.0, 900.)},
                                            bgo_kwargs={'energy_range':(350.0, 38000.0)})
        elo = (318., 7., 8.)
        ehi = (40058., 924., 913.)
        [self.assertAlmostEqual(specs[i].range[0], elo[i], places=0) for i in range(3)]
        [self.assertAlmostEqual(specs[i].range[1], ehi[i], places=0) for i in range(3)]
        

if __name__ == '__main__':
    unittest.main()