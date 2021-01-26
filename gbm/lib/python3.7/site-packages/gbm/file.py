# file.py: Module containing GBM filenaming convention and operations
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

import copy
import datetime
import os.path
import re

from .detectors import Detector
from .time import Met


class GbmFile:
    """Parse or construct a GBM standardized filename.
    
    Attributes:
        data_type (str): The datatype of the file
        detector (str): The detector with which the file is associated
        directory (str): The directory hosting the file
        extension (str): The filename extension
        meta (str): Additional metadata in the filename
        trigger (bool): True if the file is from a trigger. False otherwise
        uid (str): The unique id of the file
        version (int): The version number of the file      
    """
    REGEX_PATTERN = r'^glg_(?P<data_type>.+)_(?P<detector>[bn][0-9ab]|all)_(?P<trigger>(?:bn)?)(?P<uid>(?:\d{9}|\d{6}' \
                    r'_\d\dz|\d{6}))(?P<meta>(?:_.+)?)_v(?P<version>\d\d)\.(?P<extension>.+)$'

    def __init__(self):
        self.directory = ''
        self.trigger = False
        self.data_type = None
        self._detector = None
        self.uid = None
        self.meta = None
        self.version = 0
        self.extension = 'fit'

    def _init_by_dict(self, values):
        for key, val in values.items():
            # handle properties differently
            try:
                p = getattr(self, key)
                if isinstance(p, property):
                    p.__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute.".format(key))

    @property
    def detector(self):
        if not self._detector:
            return 'all'
        return self._detector

    @detector.setter
    def detector(self, value):
        if value == 'all':
            self._detector = None
        elif isinstance(value, Detector):
            self._detector = value
        else:
            if isinstance(value, str):
                d = Detector.from_str(value)
                self._detector = d if d else value
            elif isinstance(value, int):
                d = Detector.from_num(value)
                if d:
                    self._detector = d
                else:
                    raise ValueError("Invalid detector value")

    def version_str(self):
        """Return the file version number as a string
    
        Returns:
            str: The file version
        """
        if isinstance(self.version, int):
            v = "{:02d}".format(self.version)
        else:
            v = self.version
        return v

    def basename(self):
        """The file basename
    
        Returns:
            str: The file basename
        """
        if self.trigger:
            u = 'bn' + self.uid
        else:
            u = self.uid

        if self.meta:
            return str.format("glg_{}_{}_{}{}_v{}.{}",
                              self.data_type, self.detector, u, self.meta,
                              self.version_str(), self.extension)

        return str.format("glg_{}_{}_{}_v{}.{}",
                          self.data_type, self.detector, u, self.version_str(),
                          self.extension)

    def path(self):
        """The file path
    
        Returns:
            str: The path
        """
        return os.path.join(self.directory, self.basename())

    def __str__(self):
        return self.path()

    def __repr__(self):
        return self.path()

    @classmethod
    def create(cls, **kwargs):
        """Create a GbmFile from keywords
    
        Args:
            **kwargs: The properties of a GbmFile
        
        Returns:
            :class:`GbmFile`: The new filename object
        """
        obj = cls()
        obj._init_by_dict(kwargs)
        return obj

    @classmethod
    def from_path(cls, path):
        """Create a GbmFile from parsing a filename
    
        Args:
            path (str): A filename path
        
        Returns:
            :class:`GbmFile`: The new filename object
        """
        m = re.match(cls.REGEX_PATTERN, os.path.basename(path), re.I | re.S)

        result = None
        if m:
            result = cls.create(**m.groupdict())
            result.directory = os.path.dirname(path)

        return result

    def detector_list(self):
        """Generate a list of GbmFile objects, one for each GBM detector
    
        Returns:
            list of :class:`GbmFile`: The new filename objects
        """
        result = []
        for d in Detector:
            x = copy.copy(self)
            x.detector = d
            result.append(x)
        return result

    @classmethod
    def list_from_paths(cls, path_list, unknown=None):
        """Create a many GbmFiles from a list of filepaths
    
        Args:
            path_list (list of str): List of filepaths
        
        Returns:
            list of :class:`GbmFile`: The new filename object(s)
        """
        result = []
        for p in path_list:
            f = GbmFile.from_path(p)
            if f:
                result.append(f)
            else:
                if unknown is not None:
                    unknown.append(p)
                else:
                    raise ValueError('Unrecognized file name')
        return result


def scan_dir(path, hidden=False, recursive=False, absolute=False, regex=None):
    """
    Scans the given directory for files.

    Args:
        path (str): The root directory to scan.
        hidden (bool, optional): Set True if you want to include hidden files.
        recursive (bool, optional): Set True if you want to scan subdirectories 
                                    within the given path.
        absolute (bool, optional): Set true if you want the absolute path of 
                                   each file returned.
        regex (str): Set if you want to only return files matching the given 
                     regular expression.
    
    Yields:
        str: Full path to a file for each iteration.
    """
    for f in os.listdir(path):
        if not hidden:
            if f.startswith('.'):
                continue
        file_path = os.path.join(path, f)
        if absolute:
            file_path = os.path.abspath(file_path)
        if os.path.isfile(file_path):
            if regex and re.search(regex, f) is None:
                continue
            yield file_path
        elif recursive:
            yield from scan_dir(file_path, hidden, recursive, absolute, regex)


def all_exists(file_list, parent_dir=None):
    """
    Do all the files in the list exist in the filesystem?

    Args:
        file_list (list of str): List of file names to check
        parent_dir (str, optional): parent directory
    
    Returns:
        bool: True if all files exist
    """
    if not file_list:
        return False
    for f in file_list:
        if parent_dir is not None:
            path = os.path.join(parent_dir, f.basename())
        else:
            path = str(f)
        if not os.path.exists(path):
            return False
    return True


def has_detector(file_list, detector):
    """
    Does the file list contain a file for the given detector?

    Args:
        file_list (list of str): List of file names
        detector (str): Detector being searched
    
    Returns:
        bool: True if the list of file names includes the given detector
    """
    for f in file_list:
        if f.detector == detector:
            return True
    return False


def is_complete(file_list):
    """
    Does the file list contain a file for every detector?

    Args:
        file_list (list of str): List of files that represent a detector set
    
    Returns:
        bool: True if the file list contains a file for every detector
    """
    for d in Detector:
        if not has_detector(file_list, d):
            return False
    return True


def max_version(file_list):
    """
    Returns the maximum _version of file name in the given list

    Args:
        file_list (list of str): list of file names
    
    Returns:
        int: Largest _version number in the list
    """
    result = None
    for f in file_list:
        try:
            v = int(f.version)
            if result is None or v > result:
                result = v
        except ValueError:
            pass

    return result


def min_version(file_list):
    """
    Returns the minimum _version of file name in the given list

    Args:
        file_list (list of str): list of file names
    
    Returns:
        int: Smallest _version number in the list
    """
    result = None
    for f in file_list:
        try:
            v = int(f.version)
            if result is None or v < result:
                result = v
        except ValueError:
            pass

    return result


def ymd_path(base, name):
    if isinstance(name, str) or isinstance(name, GbmFile):
        v = name if isinstance(name, str) else name.basename()
        m = re.match(r'.*_(?:(?:bn)?)(\d{6})(?:(\d{3})|(_\d\d)?)_.*', v,
                     re.I | re.S)
        if m:
            d = datetime.datetime.strptime(m.group(1), "%y%m%d")
            return os.path.join(base, d.strftime('%Y-%m-%d'),
                                os.path.basename(name))
    elif isinstance(name, Met):
        return os.path.join(base, name.datetime.strftime('%Y-%m-%d'))
    elif isinstance(name, datetime.datetime) or isinstance(name,
                                                           datetime.date):
        return os.path.join(base, name.strftime('%Y-%m-%d'))
    raise ValueError("Can't parse a YMD value")
