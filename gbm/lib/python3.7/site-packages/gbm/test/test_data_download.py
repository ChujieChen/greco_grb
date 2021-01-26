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
import unittest, os, shutil
from gbm.finder import TriggerFtp, ContinuousFtp, TriggerCatalog, BurstCatalog

download_dir = data_dir = os.path.dirname(os.path.abspath(__file__))

class TestTriggerFtp(unittest.TestCase):
    finder = TriggerFtp()
    
    def test_set_trigger(self):
        self.finder.set_trigger('080916009')
        self.assertEqual(self.finder.num_files, 109)
        self.finder.set_trigger('170817529')
        self.assertEqual(self.finder.num_files, 128)
    
    def test_ls(self):
        self.finder.set_trigger('170817529')
        [self.assertTrue('ctime' in file) for file in self.finder.ls_ctime()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_cspec()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_rsp(ctime=False)]
        [self.assertTrue('ctime' in file) for file in self.finder.ls_rsp(cspec=False)]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_rsp2(ctime=False)]
        [self.assertTrue('ctime' in file) for file in self.finder.ls_rsp2(cspec=False)]
        [self.assertTrue('lc' in file) for file in self.finder.ls_lightcurve()]
        [self.assertTrue('.fit' in file) for file in self.finder.ls_cat_files()]
        self.assertTrue('trigdat' in self.finder.ls_trigdat()[0])
        [self.assertTrue(('healpix' in file) or ('skymap' in file) or 
                         ('loclist' in file) or ('locprob' in file) or
                         ('locplot' in file)) for file in self.finder.ls_localization()]
            
    def test_get(self):
        self.finder.set_trigger('170817529')
        self.finder.get_cat_files(download_dir)
        cat_files = self.finder.ls_cat_files()
        [os.remove(os.path.join(download_dir, file)) for file in cat_files]


class TestContinuousFtp(unittest.TestCase):
    finder = ContinuousFtp()
    
    def test_set_time(self):
        self.finder.set_time(met=604741251.0)
        self.assertEqual(self.finder.num_files, 379)
        self.finder.set_time(utc='2019-01-14T20:57:02.63')
        self.assertEqual(self.finder.num_files, 379)
        self.finder.set_time(gps=1263097406.735840)
        self.assertEqual(self.finder.num_files, 379)
    
    def test_ls(self):
        self.finder.set_time(met=604741251.0)
        [self.assertTrue('ctime' in file) for file in self.finder.ls_ctime()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_cspec()]
        [self.assertTrue('poshist' in file) for file in self.finder.ls_poshist()]
        [self.assertTrue('spechist' in file) for file in self.finder.ls_spechist()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte(full_day=True)]

    def test_get(self):
        self.finder.set_time(met=604741251.0)
        self.finder.get_poshist(download_dir)
        self.finder.get_ctime(download_dir, dets=('n0', 'n1', 'n2'))
        
        files = self.finder.ls_poshist()
        files.extend(self.finder.ls_ctime())
        for file in files:
            try:
                os.remove(os.path.join(download_dir, file))
            except:
                pass
    
    def test_reconnect(self):
        finder = ContinuousFtp()
        finder = ContinuousFtp()
        self.finder.set_time(met=604741251.0)


class TestTriggerCatalog(unittest.TestCase):
    catalog = TriggerCatalog()
    def test_attributes(self):
        self.assertEqual(self.catalog.num_cols, len(self.catalog.columns))
    
    def test_get_table(self):
        table = self.catalog.get_table()
        self.assertEqual(len(table.dtype), self.catalog.num_cols)
        table = self.catalog.get_table(columns=('trigger_name', 'ra', 'dec'))
        self.assertEqual(len(table.dtype), 3)
    
    def test_column_range(self):
        lo, hi = self.catalog.column_range('trigger_type')
        self.assertEqual(lo, 'DISTPAR')
        self.assertEqual(hi, 'UNRELOC')
        lo, hi = self.catalog.column_range('error_radius')
        self.assertEqual(lo, 0)
        self.assertEqual(hi, 93.54)
    
    def test_slice(self):
        sliced = self.catalog.slice('trigger_type', lo='GRB', hi='GRB')
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', lo=50.0, hi=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', hi=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', lo=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        
        sliced2 = self.catalog.slices([('trigger_type', 'GRB', 'GRB'),
                                       ('ra', None, 100.0), 
                                       ('dec', 20.0, None)])
        self.assertTrue(sliced2.num_rows < self.catalog.num_rows)
                     
class TestBurstCatalog(unittest.TestCase):
    catalog = BurstCatalog()
    def test_attributes(self):
        self.assertEqual(self.catalog.num_cols, len(self.catalog.columns))
    
    def test_get_table(self):
        table = self.catalog.get_table()
        self.assertEqual(len(table.dtype), self.catalog.num_cols)
        table = self.catalog.get_table(columns=('name', 'ra', 'dec'))
        self.assertEqual(len(table.dtype), 3)
    
    def test_column_range(self):
        lo, hi = self.catalog.column_range('flnc_best_fitting_model')
        self.assertEqual(lo, 'flnc_band')
        self.assertEqual(hi, 'nan')
        lo, hi = self.catalog.column_range('error_radius')
        self.assertEqual(lo, 0)
    
    def test_slice(self):
        sliced = self.catalog.slice('flnc_best_fitting_model', lo='flnc_band', 
                                    hi='flnc_band')
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', lo=50.0, hi=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', hi=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        sliced = self.catalog.slice('ra', lo=100.0)
        self.assertTrue(sliced.num_rows < self.catalog.num_rows)
        
        sliced2 = self.catalog.slices([('flnc_best_fitting_model', 'flnc_band', 'flnc_band'),
                                       ('ra', None, 100.0), 
                                       ('dec', 20.0, None)])
        self.assertTrue(sliced2.num_rows < self.catalog.num_rows)

if __name__ == '__main__':
    unittest.main()
