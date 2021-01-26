# detectors.py: Module containing the GBM detector definitions
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
from enum import Enum


class Detector(Enum):
    """The GBM detector names and pointings.
    
    Attributes:
        pointing (float, float): The spacecraft Azimuth/Zenith pointing 
        short_name (str): The short name of the detector (e.g. NaI 0 -> 'n0')    
    """
    N0 = ('NAI_00', 0, 45.89, 20.58)
    N1 = ('NAI_01', 1, 45.11, 45.31)
    N2 = ('NAI_02', 2, 58.44, 90.21)
    N3 = ('NAI_03', 3, 314.87, 45.24)
    N4 = ('NAI_04', 4, 303.15, 90.27)
    N5 = ('NAI_05', 5, 3.35, 89.79)
    N6 = ('NAI_06', 6, 224.93, 20.43)
    N7 = ('NAI_07', 7, 224.62, 46.18)
    N8 = ('NAI_08', 8, 236.61, 89.97)
    N9 = ('NAI_09', 9, 135.19, 45.55)
    NA = ('NAI_10', 10, 123.73, 90.42)
    NB = ('NAI_11', 11, 183.74, 90.32)
    B0 = ('BGO_00', 12, 0.00, 90.00)
    B1 = ('BGO_01', 13, 180.00, 90.00)

    def __init__(self, long_name, number, azimuth, zenith):
        self.long_name = long_name
        self.number = number
        self.azimuth = azimuth
        self.zenith = zenith

    def __repr__(self):
        return "Detector(\"{}\", \"{}\", {})".format(self.name, self.long_name,
                                                     self.number)

    def __str__(self):
        return str.lower(self.name)

    @property
    def short_name(self):
        return self.__str__()

    @property
    def pointing(self):
        return self.azimuth, self.zenith

    def is_nai(self):
        """Check if detector is an NaI
    
        Returns:
            bool: True if detector is NaI, False otherwise.
        """
        return self.name[0] == 'N'

    def is_bgo(self):
        """Check if detector is a BGO
    
        Returns:
            bool: True if detector is BGO, False otherwise.
        """
        return self.name[0] == 'B'

    @classmethod
    def from_str(cls, value):
        """Create a Detector from a short string name (e.g. 'n0')
        
        Args:
            value (str): The short name
        
        Returns:
            :class:`Detector`: The detector enum
        """
        if value.upper() in cls.__members__:
            return cls[value.upper()]
        # TODO: Reconsider returning None
        return None

    @classmethod
    def from_num(cls, num):
        """Create a Detector from an index number
        
        Args:
            num (int): The index number
        
        Returns:
            :class:`Detector`: The detector enum
        """
        for d in cls:
            if d.number == num:
                return d
        # TODO: Reconsider returning None
        return None

    @classmethod
    def nai(cls):
        """Get all detectors that are NaIs
    
        Returns:
            list of :class:`Detector`: The NaI detectors
        """
        return [x for x in cls if x.is_nai()]

    @classmethod
    def bgo(cls):
        """Get all detectors that are BGOs
        
        Returns:
            list of :class:`Detector`: The BGO detectors
        """
        return [x for x in Detector if x.is_bgo()]
