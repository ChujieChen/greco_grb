# tcat.py: GBM Trigger Catalog (TCAT) file class
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
import astropy.io.fits as fits
from collections import OrderedDict
from .data import DataFile


class Tcat(DataFile):
    """Class for Trigger Catalog (TCAT) files
    
    Attributes:
        datatype (str): The datatype of the file
        detector (str): The GBM detector the file is associated with
        directory (str): The directory the file is located in
        fermi_location (float, float): Fermi's orbital longitude and latitude
        filename (str): The filename
        full_path (str): The full path+filename
        headers (dict): The headers for each extension of the file
        id (str): The GBM file ID
        is_gbm_file (bool): True if the file is a valid GBM standard file, 
                            False if it is not.
        is_trigger (bool): True if the file is a GBM trigger file, False if not
        localizing_instrument (str): The localizing instrument
        location (float, float, float): RA, Dec, and localization uncertainty
        location_fermi_frame (float, float): Location in Fermi azimuth and zenith
        name (str): Name of the trigger
        time_range (float, float): The time range
        trigtime (float): The trigger time
    """

    def __init__(self):
        super(Tcat, self).__init__()
        self._headers = OrderedDict()

    @property
    def headers(self):
        return self._headers

    @property
    def time_range(self):
        return (self.headers['PRIMARY']['TSTART'],
                self.headers['PRIMARY']['TSTOP'])

    @property
    def trigtime(self):
        return self.headers['PRIMARY']['TRIGTIME']

    @property
    def location(self):
        return (self.headers['PRIMARY']['RA_OBJ'],
                self.headers['PRIMARY']['DEC_OBJ'],
                self.headers['PRIMARY']['ERR_RAD'])

    @property
    def name(self):
        return self.headers['PRIMARY']['OBJECT']

    @property
    def location_fermi_frame(self):
        return (self.headers['PRIMARY']['PHI'],
                self.headers['PRIMARY']['THETA'])

    @property
    def localizing_instrument(self):
        return self.headers['PRIMARY']['LOC_SRC']

    @property
    def fermi_location(self):
        return (self.headers['PRIMARY']['GEO_LONG'],
                self.headers['PRIMARY']['GEO_LAT'])

    @classmethod
    def open(cls, filename):
        """Open a TCAT file and return the Tcat object
        
        Args:
            filename (str):  The filename of the TCAT file
        
        Returns:
            :class:`Tcat`: The Tcat object
        """
        obj = cls()
        obj._file_properties(filename)

        # open FITS file
        with fits.open(filename) as hdulist:
            for hdu in hdulist:
                obj._headers.update({hdu.name: hdu.header})

        return obj
