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

from gbm.data.tcat import Tcat

data_dir = 'data'

class TestTcat(TestCase):
    filename = os.path.join(data_dir, 'glg_tcat_all_bn190222537_v01.fit')
    time_range = (572532678.644472, 572533293.055776)
    trigtime = 572532812.150778
    location = (147.320, 60.9400, 1.50000)
    azzen = (207.094, 142.186)
    name = 'GRB190222537'
    localizer = 'Fermi, GBM'
    lonlat = (74.2500, -24.8500)
    
    def test_attributes(self):
        t = Tcat.open(self.filename)
        self.assertEqual(len(t.headers), 1)
        self.assertEqual(t.time_range, self.time_range)
        self.assertEqual(t.trigtime, self.trigtime)
        self.assertEqual(t.location, self.location)
        self.assertEqual(t.location_fermi_frame, self.azzen)
        self.assertEqual(t.name, self.name)
        self.assertEqual(t.localizing_instrument, self.localizer)
        self.assertEqual(t.fermi_location, self.lonlat)
