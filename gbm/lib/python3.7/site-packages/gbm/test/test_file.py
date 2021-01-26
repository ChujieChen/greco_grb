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
import os
import glob

from shutil import rmtree
from subprocess import call

from gbm.file import *
from gbm.time import Met


class TestFilename(unittest.TestCase):
    @classmethod
    def create_directory_structure(cls):
        files = [

            # Sample dir structure

            "aaa/file_1.txt",
            "aaa/file_2.txt",
            "aaa/file_3.txt",
            "aaa/file_4.txt",
            "aaa/.hidden_1.txt",
            "bbb/file_1.txt",
            "bbb/file_2.txt",
            "bbb/file_3.txt",
            "bbb/file_4.txt",
            "bbb/.hidden_1.txt",
            "ccc/file_1.txt",
            "ccc/file_2.txt",
            "ccc/file_3.txt",
            "ccc/file_4.txt",
            "ccc/.hidden_1.txt",
            ".ddd/file_1.txt",
            ".ddd/file_2.txt",
            ".ddd/file_3.txt",
            ".ddd/file_4.txt",
            "complete/glg_poshist_all_160115_12z_v02.fit",
            "complete/glg_tte_n0_160115_12z_v00.fit",
            "complete/glg_tte_n1_160115_12z_v00.fit",
            "complete/glg_tte_n2_160115_12z_v00.fit",
            "complete/glg_tte_n3_160115_12z_v00.fit",
            "complete/glg_tte_n4_160115_12z_v00.fit",
            "complete/glg_tte_n5_160115_12z_v00.fit",
            "complete/glg_tte_n6_160115_12z_v00.fit",
            "complete/glg_tte_n7_160115_12z_v00.fit",
            "complete/glg_tte_n8_160115_12z_v00.fit",
            "complete/glg_tte_n9_160115_12z_v00.fit",
            "complete/glg_tte_na_160115_12z_v00.fit",
            "complete/glg_tte_nb_160115_12z_v00.fit",
            "complete/glg_tte_b0_160115_12z_v00.fit",
            "complete/glg_tte_b1_160115_12z_v00.fit",
            "incomplete/glg_tte_n0_160115_12z_v00.fit",
            "incomplete/glg_tte_n1_160115_12z_v00.fit",
            "incomplete/glg_tte_n2_160115_12z_v00.fit",
            "incomplete/glg_tte_n3_160115_12z_v00.fit",
            "incomplete/glg_tte_n4_160115_12z_v00.fit",
            "incomplete/glg_tte_n5_160115_12z_v00.fit",
            "incomplete/glg_tte_n6_160115_12z_v00.fit",
            "incomplete/glg_tte_n7_160115_12z_v00.fit",
            "incomplete/glg_tte_n9_160115_12z_v00.fit",
            "incomplete/glg_tte_na_160115_12z_v00.fit",
            "incomplete/glg_tte_nb_160115_12z_v00.fit",
            "incomplete/glg_tte_b0_160115_12z_v00.fit",
            "incomplete/glg_tte_b1_160115_12z_v00.fit",
            "_version/glg_tte_n0_160115_12z_v05.fit",
            "_version/glg_tte_n1_160115_12z_v06.fit",
            "_version/glg_tte_n2_160115_12z_v07.fit",
            "_version/glg_tte_n3_160115_12z_v08.fit",
            "_version/glg_tte_n4_160115_12z_v09.fit",
            "_version/glg_tte_n5_160115_12z_v10.fit",
            "_version/glg_tte_n6_160115_12z_v11.fit",
            "_version/glg_tte_n7_160115_12z_v12.fit",
            "_version/glg_tte_n8_160115_12z_v13.fit",
            "_version/glg_tte_n9_160115_12z_v14.fit",
            "_version/glg_tte_na_160115_12z_v15.fit",
            "_version/glg_tte_nb_160115_12z_v16.fit",
            "_version/glg_tte_b0_160115_12z_v17.fit",
            "_version/glg_tte_b1_160115_12z_v18.fit",
        ]

        def create_if_not_exists(path):
            if not os.path.isdir(path):
                os.makedirs(path)

        old_cwd = os.getcwd()
        create_if_not_exists("test_files")
        os.chdir("test_files")

        for f in files:
            directory = os.path.dirname(f)
            create_if_not_exists(directory)
            call(["touch", f])

        os.chdir(old_cwd)

    @classmethod
    def setUpClass(cls):
        cls.create_directory_structure()

    @classmethod
    def tearDownClass(cls):
        rmtree("test_files")

    def test_daily(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).ymd, data_type='test',
                           detector='n0')
        self.assertEqual(f.basename(), 'glg_test_n0_150911_v00.fit')

    def test_daily_all(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).ymd, data_type='test')
        self.assertEqual(f.basename(), 'glg_test_all_150911_v00.fit')

    def test_hourly(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).ymd_h, data_type='test',
                           detector='b1')
        self.assertEqual(f.basename(), 'glg_test_b1_150911_12z_v00.fit')

    def test_hourly_all(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).ymd_h, data_type='test')
        self.assertEqual(f.basename(), 'glg_test_all_150911_12z_v00.fit')

    def test_trig(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).bn, data_type='test',
                           detector='n0',
                           trigger=True)
        self.assertEqual(f.basename(), 'glg_test_n0_bn150911500_v00.fit')

    def test_trig_all(self):
        f = GbmFile.create(uid=Met.from_datetime(datetime.datetime(2015, 9, 11, 12)).bn, data_type='test', trigger=True)
        self.assertEqual(f.basename(), 'glg_test_all_bn150911500_v00.fit')

    def test_detector_list(self):

        expected = [
            'glg_test_n0_170514999_v00.fit',
            'glg_test_n1_170514999_v00.fit',
            'glg_test_n2_170514999_v00.fit',
            'glg_test_n3_170514999_v00.fit',
            'glg_test_n4_170514999_v00.fit',
            'glg_test_n5_170514999_v00.fit',
            'glg_test_n6_170514999_v00.fit',
            'glg_test_n7_170514999_v00.fit',
            'glg_test_n8_170514999_v00.fit',
            'glg_test_n9_170514999_v00.fit',
            'glg_test_na_170514999_v00.fit',
            'glg_test_nb_170514999_v00.fit',
            'glg_test_b0_170514999_v00.fit',
            'glg_test_b1_170514999_v00.fit'

        ]

        f = GbmFile.create(uid='170514999', version=0, data_type='test')

        file_list = f.detector_list()

        result = [str(x) for x in file_list]
        self.assertEqual(expected, result)

    def test_from_path_triggered(self):
        f = "fake/glg_tte_n5_bn171115500_v05.fit"
        g = GbmFile.from_path(f)
        self.assertIsNotNone(g)
        self.assertTrue(g.trigger)
        self.assertEqual(f, g.path())

    def test_from_path_continuous(self):
        f = "fake/glg_tte_n5_171115500_v05.fit"
        g = GbmFile.from_path(f)
        self.assertIsNotNone(g)
        self.assertFalse(g.trigger)
        self.assertEqual(f, g.path())

    def test_from_path_scat(self):
        f = "fake/glg_scat_all_bn080714086_flnc_band_v00.fit"
        g = GbmFile.from_path(f)
        self.assertIsNotNone(g)
        self.assertEqual('scat', g.data_type)
        self.assertEqual('all', g.detector)
        self.assertTrue(g.trigger)
        self.assertEqual('080714086', g.uid)
        self.assertEqual('_flnc_band', g.meta)
        self.assertEqual('00', g.version)
        self.assertEqual('fit', g.extension)
        self.assertEqual(f, g.path())

    def test_scan_dir(self):
        expected = [
            "test_files/aaa/file_1.txt",
            "test_files/aaa/file_2.txt",
            "test_files/aaa/file_3.txt",
            "test_files/aaa/file_4.txt",
            "test_files/bbb/file_1.txt",
            "test_files/bbb/file_2.txt",
            "test_files/bbb/file_3.txt",
            "test_files/bbb/file_4.txt",
            "test_files/ccc/file_1.txt",
            "test_files/ccc/file_2.txt",
            "test_files/ccc/file_3.txt",
            "test_files/ccc/file_4.txt",
            "test_files/complete/glg_poshist_all_160115_12z_v02.fit",
            "test_files/complete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n8_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/complete/glg_tte_na_160115_12z_v00.fit",
            "test_files/complete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_na_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/_version/glg_tte_n0_160115_12z_v05.fit",
            "test_files/_version/glg_tte_n1_160115_12z_v06.fit",
            "test_files/_version/glg_tte_n2_160115_12z_v07.fit",
            "test_files/_version/glg_tte_n3_160115_12z_v08.fit",
            "test_files/_version/glg_tte_n4_160115_12z_v09.fit",
            "test_files/_version/glg_tte_n5_160115_12z_v10.fit",
            "test_files/_version/glg_tte_n6_160115_12z_v11.fit",
            "test_files/_version/glg_tte_n7_160115_12z_v12.fit",
            "test_files/_version/glg_tte_n8_160115_12z_v13.fit",
            "test_files/_version/glg_tte_n9_160115_12z_v14.fit",
            "test_files/_version/glg_tte_na_160115_12z_v15.fit",
            "test_files/_version/glg_tte_nb_160115_12z_v16.fit",
            "test_files/_version/glg_tte_b0_160115_12z_v17.fit",
            "test_files/_version/glg_tte_b1_160115_12z_v18.fit",

        ]

        files = [x for x in scan_dir("test_files", recursive=True)]

        files.sort()
        expected.sort()

        self.assertListEqual(expected, files)

    def test_scan_dir_absolute(self):
        expected = [
            "test_files/aaa/file_1.txt",
            "test_files/aaa/file_2.txt",
            "test_files/aaa/file_3.txt",
            "test_files/aaa/file_4.txt",
            "test_files/bbb/file_1.txt",
            "test_files/bbb/file_2.txt",
            "test_files/bbb/file_3.txt",
            "test_files/bbb/file_4.txt",
            "test_files/ccc/file_1.txt",
            "test_files/ccc/file_2.txt",
            "test_files/ccc/file_3.txt",
            "test_files/ccc/file_4.txt",
            "test_files/complete/glg_poshist_all_160115_12z_v02.fit",
            "test_files/complete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n8_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/complete/glg_tte_na_160115_12z_v00.fit",
            "test_files/complete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_na_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/_version/glg_tte_n0_160115_12z_v05.fit",
            "test_files/_version/glg_tte_n1_160115_12z_v06.fit",
            "test_files/_version/glg_tte_n2_160115_12z_v07.fit",
            "test_files/_version/glg_tte_n3_160115_12z_v08.fit",
            "test_files/_version/glg_tte_n4_160115_12z_v09.fit",
            "test_files/_version/glg_tte_n5_160115_12z_v10.fit",
            "test_files/_version/glg_tte_n6_160115_12z_v11.fit",
            "test_files/_version/glg_tte_n7_160115_12z_v12.fit",
            "test_files/_version/glg_tte_n8_160115_12z_v13.fit",
            "test_files/_version/glg_tte_n9_160115_12z_v14.fit",
            "test_files/_version/glg_tte_na_160115_12z_v15.fit",
            "test_files/_version/glg_tte_nb_160115_12z_v16.fit",
            "test_files/_version/glg_tte_b0_160115_12z_v17.fit",
            "test_files/_version/glg_tte_b1_160115_12z_v18.fit",

        ]

        files = [x for x in scan_dir("test_files", absolute=True, recursive=True)]

        files.sort()
        expected.sort()
        for f in files:
            e = os.path.abspath(expected.pop(0))
            self.assertEqual(f, e)

    def test_scan_dir_hidden(self):
        expected = [
            "test_files/aaa/file_1.txt",
            "test_files/aaa/file_2.txt",
            "test_files/aaa/file_3.txt",
            "test_files/aaa/file_4.txt",
            "test_files/aaa/.hidden_1.txt",
            "test_files/bbb/file_1.txt",
            "test_files/bbb/file_2.txt",
            "test_files/bbb/file_3.txt",
            "test_files/bbb/file_4.txt",
            "test_files/bbb/.hidden_1.txt",
            "test_files/ccc/file_1.txt",
            "test_files/ccc/file_2.txt",
            "test_files/ccc/file_3.txt",
            "test_files/ccc/file_4.txt",
            "test_files/ccc/.hidden_1.txt",
            "test_files/.ddd/file_1.txt",
            "test_files/.ddd/file_2.txt",
            "test_files/.ddd/file_3.txt",
            "test_files/.ddd/file_4.txt",
            "test_files/complete/glg_poshist_all_160115_12z_v02.fit",
            "test_files/complete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n8_160115_12z_v00.fit",
            "test_files/complete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/complete/glg_tte_na_160115_12z_v00.fit",
            "test_files/complete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/complete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n1_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n2_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n3_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n4_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n5_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n6_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n7_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_n9_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_na_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_nb_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b0_160115_12z_v00.fit",
            "test_files/incomplete/glg_tte_b1_160115_12z_v00.fit",
            "test_files/_version/glg_tte_n0_160115_12z_v05.fit",
            "test_files/_version/glg_tte_n1_160115_12z_v06.fit",
            "test_files/_version/glg_tte_n2_160115_12z_v07.fit",
            "test_files/_version/glg_tte_n3_160115_12z_v08.fit",
            "test_files/_version/glg_tte_n4_160115_12z_v09.fit",
            "test_files/_version/glg_tte_n5_160115_12z_v10.fit",
            "test_files/_version/glg_tte_n6_160115_12z_v11.fit",
            "test_files/_version/glg_tte_n7_160115_12z_v12.fit",
            "test_files/_version/glg_tte_n8_160115_12z_v13.fit",
            "test_files/_version/glg_tte_n9_160115_12z_v14.fit",
            "test_files/_version/glg_tte_na_160115_12z_v15.fit",
            "test_files/_version/glg_tte_nb_160115_12z_v16.fit",
            "test_files/_version/glg_tte_b0_160115_12z_v17.fit",
            "test_files/_version/glg_tte_b1_160115_12z_v18.fit",

        ]

        files = [x for x in scan_dir("test_files", recursive=True, hidden=True)]

        files.sort()
        expected.sort()

        self.assertListEqual(expected, files)

    def test_hidden_regex(self):
        expected = [
            "test_files/aaa/.hidden_1.txt",
            "test_files/bbb/.hidden_1.txt",
            "test_files/ccc/.hidden_1.txt",
        ]

        files = [x for x in scan_dir("test_files", recursive=True, hidden=True, regex="hidden")]

        files.sort()
        expected.sort()

        self.assertListEqual(expected, files)

    def test_is_complete(self):
        dt = datetime.datetime(2016, 1, 15, 12, 0)
        glob_expression = GbmFile.create(data_type='tte', uid=Met.from_datetime(dt).ymd_h, detector='*').basename()
        files = glob.glob(os.path.join('test_files', 'complete', glob_expression))
        fn_list = GbmFile.list_from_paths(files)
        self.assertEqual(is_complete(fn_list), True)
        files = glob.glob(os.path.join('test_files', 'incomplete', glob_expression))
        fn_list = GbmFile.list_from_paths(files)
        self.assertEqual(is_complete(fn_list), False)

    def test_max_min_versions(self):
        dt = datetime.datetime(2016, 1, 15, 12, 0)
        glob_expression = GbmFile.create(data_type='tte', uid=Met.from_datetime(dt).ymd_h, detector='*',
                                         version='??').basename()
        files = glob.glob(os.path.join('test_files', '_version', glob_expression))
        fn_list = GbmFile.list_from_paths(files)
        self.assertEqual(max_version(fn_list), 18)
        self.assertEqual(min_version(fn_list), 5)

    def test_all_exists(self):
        dt = datetime.datetime(2016, 1, 15, 12, 0)
        files = GbmFile.create(data_type='tte', uid=Met.from_datetime(dt).ymd_h, detector=Detector.N0,
                               version=0).detector_list()
        self.assertEqual(all_exists(files, os.path.join('test_files', 'complete')), True)
        self.assertEqual(all_exists(files, os.path.join('test_files', 'incomplete')), False)

    def test_ymd_path_hourly_str(self):
        b = '/data'
        f = 'glg_tte_n5_160712_12z_v00.fit'
        p = ymd_path(b, f)
        self.assertEqual(p, '/data/2016-07-12/glg_tte_n5_160712_12z_v00.fit')

    def test_ymd_path_trig_str(self):
        b = '/data'
        f = 'glg_tte_n0_bn170710500_v00.fit'
        p = ymd_path(b, f)
        self.assertEqual(p, '/data/2017-07-10/glg_tte_n0_bn170710500_v00.fit')

    def test_ymd_path_daily_str(self):
        b = '/data'
        f = 'glg_poshist_all_170802_v00.fit'
        p = ymd_path(b, f)
        self.assertEqual(p, '/data/2017-08-02/glg_poshist_all_170802_v00.fit')

    def test_ymd_path_legacy_str(self):
        b = '/data'
        f = 'glg_ctime_na_150413125_v00.pha'
        p = ymd_path(b, f)
        self.assertEqual(p, '/data/2015-04-13/glg_ctime_na_150413125_v00.pha')

if __name__ == '__main__':
    unittest.main()