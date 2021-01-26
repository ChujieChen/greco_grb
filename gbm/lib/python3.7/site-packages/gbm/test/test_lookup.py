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
from gbm.lookup.lookup import *


class TestLookup(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestLookup, self).__init__(*args, **kwargs)
        self.datafiles = [
            {
                "filename": "glg_cspec_b1_bn090926181_v00.pha",
                "response": None,
                "background": {
                    "method": "Polynomial",
                    "datatype": "binned",
                    "args": [2],
                    "kwargs": {}
                },
                "binnings": {
                    "time": None,
                    "energy": None,
                },
                "selections": {
                    "source": [[0.7054115, 19.058352]],
                    "energy": [[250.0, 38000.0]],
                    "background": [[69.57779188816323, 240.8593380212668], [-45.687147962415935, -4.751935679032698]],
                },
                "views": {
                    "time": {
                        'xmin': -14.447117535623448,
                        'xmax': 50.187428174981676,
                        'ymin': 1831.095164143789,
                        'ymax': 5137.280091318506
                    },
                    "energy": None
                }
            },
            {
                "filename": "glg_cspec_n6_bn090926181_v00.pha",
                "response": "glg_cspec_n6_bn090926181_v00.rsp2",
                "background": {
                    "method": "Polynomial",
                    "datatype": "binned",
                    "args": [3],
                    "kwargs": {}
                },
                "binnings": {
                    "time": None,
                    "energy": None,
                },
                "selections": {
                    "source": [[0.7054115, 19.058352]],
                    "energy": [[8.47545, 905.397]],
                    "background": [[-1529.48, -1200.07], [-53.5312, -17.6488], [141.106, 1764.64]],
                },
                "views": {
                    "time": {
                        'xmin': -10.0,
                        'xmax': 50.0,
                        'ymin': 0.0,
                        'ymax': 9556.3
                    },
                    "energy": {
                        'xmin': 4.44444,
                        'xmax': 2000.0,
                        'ymin': 0.0352183,
                        'ymax': 69.4991
                    }
                }
            }
        ]

    def assert_datafile(self, df, d):
        self.assertEqual(df.filename, d.get('filename', None))
        self.assertEqual(df.response, d.get('response', None))

        bkg = d.get('background', None)
        if bkg is None:
            self.assertIsNone(df.background)
        else:
            self.assertIsNotNone(df.background)
            self.assertEqual(df.background.method, bkg.get('method', None))
            self.assertEqual(list(df.background.args), bkg.get('args', None))
            self.assertEqual(df.background.kwargs, bkg.get('kwargs', None))
            self.assertEqual(df.background.datatype, bkg.get('datatype', None))

        if 'binnings' in d:
            b = d['binnings']
            self.assertEqual(df.binnings.energy, b.get('energy', None))
            self.assertEqual(df.binnings.time, b.get('time', None))
        else:
            self.assertIsNone(df.binnings.energy)
            self.assertIsNone(df.binnings.time)

        if 'selections' in d:
            b = d['selections']
            self.assertEqual(df.selections.background, b.get('background', None))
            self.assertEqual(df.selections.energy, b.get('energy', None))
            self.assertEqual(df.selections.source, b.get('source', None))
        else:
            self.assertIsNone(df.selections.background)
            self.assertIsNone(df.selections.energy)
            self.assertIsNone(df.selections.source)

        if 'views' in d:
            b = d['views']
            v = b.get('energy', None)
            if v is None:
                self.assertIsNone(df.views.energy)
            else:
                self.assertIsInstance(df.views.energy, View)
                self.assertEqual(df.views.energy, View.from_dict(v))
            v = b.get('time', None)
            if v is None:
                self.assertIsNone(df.views.time)
            else:
                self.assertIsInstance(df.views.time, View)
                self.assertEqual(df.views.time, View.from_dict(v))
        else:
            self.assertIsNone(df.views.energy)
            self.assertIsNone(df.views.time)

    def test_create_datafile_from_dict(self):
        for d in self.datafiles:
            df = DataFileLookup.from_dict(d)
            self.assert_datafile(df, d)

    def test_create_datafile(self):
        """Create a DataFile object with:
        {
            "filename": "glg_cspec_n6_bn090926181_v00.pha",
            "response": "glg_cspec_n6_bn090926181_v00.rsp2",
            "background": {
                "method": "Polynomial",
                "datatype": "binned",
                "args": [ 3 ],
                "kwargs": {}
            },
            "time_binning": null,
            "energy_binning": null,
            "source_selection": [[0.7054115, 19.058352]],
            "energy_selection": [[8.47545, 905.397]],
            "background_selection": [[-1529.48, -1200.07], [-53.5312, -17.6488], [141.106, 1764.64]],
            "time_display_view": [-10.0, 50.0, 0.0, 9556.3],
            "energy_display_view": [4.44444, 2000.0, 0.0352183, 69.4991]
        }
        """
        df = DataFileLookup()
        df.filename = "glg_cspec_n6_bn090926181_v00.pha"
        df.response = "glg_cspec_n6_bn090926181_v00.rsp2"

        df.background = LookupBackground()
        df.background.method = "Polynomial"
        df.background.datatype = "binned"
        df.background.args = (3, )
        df.background.kwargs = {}

        df.selections.source = [[0.7054115, 19.058352]]
        df.selections.energy = [[8.47545, 905.397]]
        df.selections.background = [[-1529.48, -1200.07], [-53.5312, -17.6488], [141.106, 1764.64]]
        df.views.time = View(-10.0, 50.0, 0.0, 9556.3)
        df.views.energy = View(4.44444, 2000.0, 0.0352183, 69.4991)

        self.assert_datafile(df, self.datafiles[1])

    def test_read_write_lookup(self):
        lu_w = LookupFile()
        for d in self.datafiles:
            df = DataFileLookup.from_dict(d)
            lu_w[df.filename] = df
        lu_w.write_to("test.json")

        lu_r = LookupFile.read_from("test.json")

        self.assertEqual(lu_r.file_date, lu_w.file_date)
        self.assertEqual(len(lu_r.datafiles), len(lu_w.datafiles))
        for k, df_r in lu_r.datafiles.items():
            df_w = lu_w.datafiles[k]

            self.assertEqual(df_r.filename, df_w.filename)
            self.assertEqual(df_r.response, df_w.response)

            self.assertEqual(df_r.binnings.energy, df_w.binnings.energy)
            self.assertEqual(df_r.binnings.time, df_w.binnings.time)

            self.assertEqual(df_r.selections.background, df_w.selections.background)
            self.assertEqual(df_r.selections.energy, df_w.selections.energy)
            self.assertEqual(df_r.selections.source, df_w.selections.source)

            self.assertEqual(df_r.views.energy, df_w.views.energy)
            self.assertEqual(df_r.views.time, df_w.views.time)

    def test_read_from_rmfit(self):
        lu = LookupFile.read_from_rmfit('data/glg_ctime_nb_bn120415958_v00.lu')
        df = lu.datafiles['glg_ctime_nb_bn120415958_v00.pha']

        # TODO: Finalize the background format for the Lookup file
        self.assertEqual(df.background.method, "Polynomial")
        self.assertEqual(df.background.args[0], 2)
        
        self.assertEqual(df.views.time, View(-60.0, 40.0, 400.0, 1400.0))
        self.assertEqual(df.views.energy, View(4.32375, 2000.0, 0.0676773, 22.8421))

        self.assertIsNone(df.binnings.energy)
        time_edges = df.binnings.time[0].args[0]
        self.assertEqual(1861, len(time_edges))
        self.assertEqual(0, time_edges[0])
        self.assertEqual(14426, time_edges[-1])

        energy_select = df.selections.energy
        self.assertEqual(len(energy_select), 1)
        self.assertEqual(18.8469, energy_select[0][0])
        self.assertEqual(767.693, energy_select[0][1])

        time_select = df.selections.source
        self.assertEqual(len(time_select), 2)
        self.assertEqual(time_select[0][0], -17.647648)
        self.assertEqual(time_select[0][1], -14.902550)
        self.assertEqual(time_select[1][0], -5.2947070)
        self.assertEqual(time_select[1][1], 2.7445084)

        bkgd_select = df.selections.background
        self.assertEqual(len(bkgd_select), 2)
        self.assertEqual(bkgd_select[0][0], -480.00706)
        self.assertEqual(bkgd_select[0][1], -115.30119)
        self.assertEqual(bkgd_select[1][0], 72.934103)
        self.assertEqual(bkgd_select[1][1], 522.34586)

    def test_read_ti(self):
        lu = LookupFile.read_from_rmfit('data/glg_tte_n9_bn090131090_v00.lu', ti_file='data/glg_tte_n9_bn090131090_v00.ti')

        tte_edges = lu.datafiles['glg_tte_n9_bn090131090_v00.fit'].binnings.time[0].args[0]
        self.assertEqual(len(tte_edges), 319)
        self.assertEqual(tte_edges[0], -25.600000)
        self.assertEqual(tte_edges[-1], 300.73552)

if __name__ == '__main__':
    unittest.main()
      
