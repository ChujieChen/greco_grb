# data.py: Data file class definition
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

import os
from gbm.file import GbmFile
from gbm.time import Met


class DataFile:
    """Base class for an interface to data files
    
    Note:
        This class should not be used directly, instead use one that inherits 
        from it and is specific to your data type e.g. :class:`~gbm.data.TTE`
    
    Attributes:
        datatype (str): The datatype of the file
        detector (str): The GBM detector the file is associated with
        directory (str): The directory the file is located in
        filename (str): The filename
        full_path (str): The full path+filename
        id (str): The GBM file ID
        is_gbm_file (bool): True if the file is a valid GBM standard file, 
                            False if it is not.
        is_trigger (bool): True if the file is a GBM trigger file, False if not
    """
    def __init__(self):
        self._full_path = None
        self._dir = None
        self._filename = None
        self._is_gbm_file = None
        self._filename_obj = None
    
    def __str__(self):
        return self.filename

    def _file_properties(self, filename, ignore_file_check=False):
        
        if not ignore_file_check:
            if not os.path.isfile(filename):
                raise IOError("File {0} does not exist".format(filename))
    
        self._full_path = filename
        self._dir = os.path.dirname(filename)
        self._filename = os.path.basename(filename)
        self._is_gbm_file = False
        try:
            self._filename_obj = GbmFile.from_path(filename)
            self._is_gbm_file = True
        except:
            pass
    
    @property
    def is_gbm_file(self):
        return self._is_gbm_file
    
    @property
    def id(self):
        if self.is_gbm_file:
            return self._filename_obj.uid
    
    @property
    def filename(self):
        return self._filename
    
    @property
    def is_trigger(self):
        if self.is_gbm_file:
            if self._filename_obj.trigger == 'bn':
                return True
        return False
    
    @property
    def detector(self):
        if self.is_gbm_file:
            try:
                return self._filename_obj.detector.short_name
            except:
                return 'all'
    
    @property
    def datatype(self):
        if self.is_gbm_file:
            return self._filename_obj.data_type.upper()
    
    @property
    def directory(self):
        return self._dir
    
    @property
    def full_path(self):
        return self._full_path
    
    def set_properties(self, detector=None, trigtime=None, tstart=None, 
                       datatype=None, extension=None, meta=None):        
        """Set the properties of the data file.  Useful for creating new
        data files.  A standardized filename will be built from the 
        parameters
        
        Args:
            detector (str, optional): The detector the data file belongs to
            trigtime (float, optional): The trigger time, if applicable
            tstart (float):  The start time of the data file. 
                             Must be set if trigtime is not.
            datatype (str, optional): The type of data the file contains
            extension (str, optional): The extension of the data file
            meta (str, optional): A metadata atttribute to be added to the filename
        """
        if (trigtime is None) and (tstart is None):
            raise KeyError('Either trigtime or tstart need to be defined')
        
        if trigtime is not None:
            trigger = True
            met = Met(trigtime)
            id = met.bn
        else:
            trigger = False
            met = Met(tstart)
            id = met.ymd_h
            
        filename = GbmFile.create(uid=id, data_type=datatype, detector=detector, 
                                  trigger=trigger, extension=extension, meta=meta)
        
        self._file_properties(filename.basename(), ignore_file_check=True)
