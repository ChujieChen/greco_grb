# lookup.py: GSpec and RMfit lookup classes
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
import datetime as dt
import json
import os.path
import warnings

import numpy as np

from gbm.detectors import Detector
from gbm.file import GbmFile
from gbm.types import ListReader


class LookupMethod:
    """Defines the attributes of a method call"""

    def __init__(self):
        self.method = None
        self.args = None
        self.kwargs = {}

    @classmethod
    def from_dict(cls, d):
        r = cls()
        r.method = d.get('method', None)
        r.args = tuple(d.get('args', None))
        r.kwargs = d.get('kwargs', {})
        return r


class LookupBackground(LookupMethod):
    """Defines the attributes of a background binning method"""

    def __init__(self):
        super(LookupBackground, self).__init__()
        self.datatype = None

    @classmethod
    def from_dict(cls, d):
        r = super(LookupBackground, cls).from_dict(d)
        r.datatype = d.get('datatype', None)
        return r


class LookupEnergyBinning(LookupMethod):
    """Defines the attributes of an energy binning method"""

    def __init__(self):
        super(LookupEnergyBinning, self).__init__()
        self.start = None
        self.stop = None

    @classmethod
    def from_dict(cls, d):
        r = super(LookupEnergyBinning, cls).from_dict(d)
        r.start = d.get('start', None)
        r.stop = d.get('stop', None)
        return r


class LookupTimeBinning(LookupEnergyBinning):
    """Defines the attributes of a time binning method"""

    def __init__(self):
        super(LookupTimeBinning, self).__init__()
        self.datatype = None

    @classmethod
    def from_dict(cls, d):
        r = super(LookupTimeBinning, cls).from_dict(d)
        r.datatype = d.get('datatype', None)
        return r


class View:
    """Defines the bounds of a view"""

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __eq__(self, other):
        return self.xmin == other.xmin and self.xmax == other.xmax and self.ymin == other.ymin \
               and self.ymax == other.ymax

    @classmethod
    def from_dict(cls, d):
        r = cls()
        if d:
            r.xmin = d.get('xmin', None)
            r.xmax = d.get('xmax', None)
            r.ymin = d.get('ymin', None)
            r.ymax = d.get('ymax', None)
        return r

    @classmethod
    def from_list(cls, l):
        r = cls()
        r.xmin = l[0]
        r.xmax = l[1]
        r.ymin = l[2]
        r.ymax = l[3]
        return r

    def to_list(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax]

    def xrange(self):
        return self.xmin, self.xmax

    def yrange(self):
        return self.ymin, self.ymax


class Binnings:

    def __init__(self):
        self.energy = None
        self.time = None


class Selections:

    def __init__(self):
        self.background = None
        self.energy = None
        self.source = None

    def add(self, type, item):
        if getattr(self, type) is None:
            setattr(self, type, list())
        getattr(self, type).append(item)


class Views:

    def __init__(self):
        self.energy = None
        self.time = None


class DataFileLookup:
    """Defines all the information associated with a datafile"""

    def __init__(self):
        self.filename = None
        self.detector = None
        self.response = None
        self.background = None
        self.binnings = Binnings()
        self.selections = Selections()
        self.views = Views()

    @classmethod
    def from_dict(cls, d):

        def set_attributes(obj, d):
            for k, v in d.items():
                setattr(obj, k, v)

        r = cls()
        r.filename = d.get('filename', None)

        det = d.get('detector', None)
        if det:
            r.detector = Detector.from_str(det)

        r.response = d.get('response', None)

        bkg = d.get('background', None)
        if bkg:
            r.background = LookupBackground.from_dict(bkg)

        binnings = d.get('binnings', None)
        if binnings:
            energies = binnings.get('energy', None)
            if energies:
                for e in energies:
                    if r.binnings.energy is None:
                        r.binnings.energy = [LookupEnergyBinning.from_dict(e)]
                    else:
                        r.binnings.energy.append(
                            LookupEnergyBinning.from_dict(e))
            times = binnings.get('time', None)
            if times:
                for t in times:
                    if r.binnings.time is None:
                        r.binnings.time = [LookupTimeBinning.from_dict(t)]
                    else:
                        r.binnings.time.append(LookupTimeBinning.from_dict(t))

        if 'selections' in d:
            set_attributes(r.selections, d['selections'])

        views = d.get('views', None)
        if views:
            e = views.get('energy', None)
            if e:
                r.views.energy = View.from_dict(e)
            t = views.get('time', None)
            if t:
                r.views.time = View.from_dict(t)
        return r

    @staticmethod
    def assert_selections(selections):
        """Check to ensure the selections are of the correct form.

        Parameters:
        -----------
        selections: tuple or list of tuples
            The selection(s) to check

        Returns:
        --------
        selections: list
        """
        if (all(isinstance(selection, list) for selection in selections)) | \
                (
                all(isinstance(selection, tuple) for selection in selections)):
            if any(len(selection) != 2 for selection in selections):
                raise ValueError('Each range in selections must be of the '
                                 'form (lo, hi)')
            else:
                return selections
        else:
            if len(selections) != 2:
                raise ValueError('Selections must either be a range of '
                                 'the form (lo, hi) or a list of ranges')
            else:
                return [selections]

    def set_response(self, rsp_filename):
        """Add a response file for the data

        Parameters:
        --------------
        rsp_filename: str
            The filename of the response file
        """
        if rsp_filename is None:
            self.response = None
        else:
            self.response = os.path.basename(rsp_filename)

    def set_background_model(self, background_name, datatype, *args, **kwargs):
        """Add a new background model for the data file

        Parameters:
        --------------
        background_class: str
            The background fitting/estimation name
        datatype: str
            The datatype the background is applied to. Either 'binned' or 'unbinned'
        *args:
            Additional arguments used by the background class
        **kwargs:
            Additional keywords used by the background class
        """
        bkg = LookupBackground()
        bkg.method = background_name
        bkg.datatype = datatype
        bkg.args = args
        bkg.kwargs = kwargs
        self.background = bkg

    def set_time_binning(self, binning_name, datatype, *args, start=None,
                         stop=None, **kwargs):
        """Add a new time binning function for the data file

        Parameters:
        --------------
        binning_function: str
            The binning function name
        datatype: str
            The datatype the binning is applied to. Either 'binned' or 'unbinned'
        *args:
            Additional arguments used by the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the
            beginning of the histogram segment.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the histogram segment.
        **kwargs:
            Additional keywords used by the binning function
        """
        time_bin = LookupTimeBinning()
        time_bin.method = binning_name
        time_bin.datatype = datatype
        time_bin.args = args
        time_bin.start = start
        time_bin.stop = stop
        time_bin.kwargs = kwargs

        if self.binnings.time is None:
            self.binnings.time = [time_bin]
        else:
            self.binnings.time.append(time_bin)

    def set_energy_binning(self, binning_function, *args, start=None,
                           stop=None, **kwargs):
        """Add a new energy binning function for the data file

        Parameters:
        --------------
        binning_function: function
            The binning function
        *args:
            Additional arguments used by the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the
            beginning of the histogram segment.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the histogram segment.
        **kwargs:
            Additional keywords used by the binning function
        """
        energy_bin = LookupEnergyBinning()
        energy_bin.method = binning_function
        energy_bin.args = args
        energy_bin.start = start
        energy_bin.stop = stop
        energy_bin.kwargs = kwargs

        if self.binnings.energy is None:
            self.binnings.energy = [energy_bin]
        else:
            self.binnings.energy.append(energy_bin)

    def set_source_selection(self, source_intervals):
        """Add source selection(s) for the data file

        Parameters:
        --------------
        dataname: str
            The data filename
        source_intervals: list
            A list of source selection intervals, each item of the list being a tuple
            of the format (low, high)
        """
        source_intervals = self.assert_selections(source_intervals)
        self.selections.source = source_intervals

    def set_energy_selection(self, energy_intervals):
        """Add energy selection(s) for the data file

        Parameters:
        --------------
        energy_intervals: list
            A list of energy selection intervals, each item of the list being a tuple
            of the format (low, high)
        """
        energy_intervals = self.assert_selections(energy_intervals)
        self.selections.energy = energy_intervals

    def set_background_selection(self, background_intervals):
        """Add background selection(s) for the data file

        Parameters:
        --------------
        background_intervals: list
            A list of background selection intervals, each item of the list being a tuple
            of the format (low, high)
        """
        self.selections.background = background_intervals

    def add_time_display_view(self, display_range):
        """Add the display range of the lightcurve for the data file

        Parameters:
        --------------
        display_range: list
            The values of the lightcurve display window in the format
            [xmin, xmax, ymin, ymax]
        """
        self.views.time = View(display_range[0], display_range[1],
                               display_range[2], display_range[3])

    def add_energy_display_view(self, display_range):
        """Add the display range of the count spectrum for the data file

        Parameters:
        --------------
        dataname: str
            The data filename
        display_range: list
            The values of the count spectrum display window in the format
            [xmin, xmax, ymin, ymax]
        """
        self.views.energy = View(display_range[0], display_range[1],
                                 display_range[2], display_range[3])


class LookupFile:
    """Class for an Gspec lookup file

    The lookup file contains one or more data files.
    """

    def __init__(self, *args, **kwargs):
        self.file_date = None
        self.datafiles = dict()

    def __getitem__(self, name):
        return self.datafiles[name]

    def __delitem__(self, key):
        del self.datafiles[key]

    def __setitem__(self, key, value):
        if isinstance(value, DataFileLookup):
            self.datafiles[key] = value
        else:
            raise ValueError("not a DataFile")

    def files(self):
        """Return the data filenames contained within the lookup"""
        return self.datafiles.keys()

    def assert_has_datafile(self, dataname):
        """Check to see if the data file has been added to the lookup

        Parameters:
        --------------
        dataname: str
            The data file name
        """
        if dataname not in self.datafiles.keys():
            raise KeyError('File {0} not currently tracked. Add this file to '
                           'the lookup and try again.'.format(dataname))

    def add_data_file(self, filepath):
        df = DataFileLookup()
        fn = GbmFile.from_path(filepath)
        df.filename = fn.basename()
        df.detector = fn.detector
        self.datafiles[df.filename] = df

    @classmethod
    def from_dict(cls, d):
        r = cls()
        r.file_date = d.get('file_date', None)
        datafiles = d.get('datafiles', None)
        if datafiles:
            for k, v in datafiles.items():
                df = DataFileLookup.from_dict(v)
                df.filename = k
                r.datafiles[k] = df
        return r

    def write_to(self, fpath):
        """
        Write contents of LookupFile to the given file path as a JSON file.
        :param fpath: full pathname for JSON file
        :return: None
        """
        self.file_date = dt.datetime.utcnow().isoformat()
        with open(fpath, "w") as fp:
            json.dump(self, fp, cls=LookupEncoder, indent=4)

    @classmethod
    def read_from(cls, fpath):
        """
        Load values to LookupFile from the JSON file at the given path.
        :param fpath: full pathname for JSON file
        :return: new LookupFile object
        """
        with open(fpath, "r") as fp:
            j = json.load(fp)
        return cls.from_dict(j)

    @classmethod
    def read_from_rmfit(cls, fpath, ti_file=None, dataname=None):
        """
        Load values to LookupFile from the RMFIT created lookup file at the given path.
        :param fpath: full pathname for RMFIT created lookup file
        :param ti_file: full pathname for RMFIT created ti file.
        :param dataname: the name of the datafile to associate this lookup file with
        :return: new LookupFile object
        """

        # RMFit selections are an 2xN array where the first element is the start values and the second element
        # are the end values. It needs to be transposed into a Nx2 array. Drop first is used to drop the convex
        # hull if the selections contain one.
        def transform_selections(x, drop_first=False):
            result = None
            if x:
                result = np.array(x).reshape(2, -1).transpose().tolist()
                if drop_first:
                    result = result[1:]
            return result

        # Begin loading RMFit lookup file making the contents a list of tokens.
        tokens = []
        with open(fpath, 'r') as contents:
            for line in contents:
                x = line.strip().split()
                if x:
                    try:
                        # If the first element a number? Then add the array to the tokens.
                        float(x[0])
                        tokens += x
                    except ValueError:
                        # Otherwise, it's a string and we will append the entire line as a token.
                        tokens.append(line)

        # Let's create the DataFile object
        data_file = DataFileLookup()

        # The input data file is based on the lookup filename
        f = GbmFile.from_path(fpath)

        if f.extension == 'lu':
            if f.data_type == 'ctime' or f.data_type == 'cspec':
                f.extension = 'pha'
            elif f.data_type == 'tte':
                f.extension = 'fit'
            else:
                raise ValueError('Not a valid lookup filename')
        else:
            raise ValueError("Not a valid lookup filename")

        if dataname:
            data_file.filename = os.path.basename(dataname)
        else:
            data_file.filename = f.basename()

        lr = ListReader(tokens)

        # energy edges, if None is returned we need to read the next value anyway which should be zero.
        energy_edges = lr.get_n(int, rmfit=True)
        if energy_edges:
            # TODO: add_energy_binning unresolved for class 'DataFileLookup'
            data_file.add_energy_binning('By Edge Index',
                                         np.array(energy_edges))

        # energy selections, if None is returned we need to read the next value anyway which should be zero.
        data_file.selections.energy = transform_selections(
            lr.get_n(float, rmfit=True), drop_first=True)

        # rebinned time edges, if None is returned we need to read the next value anyway which should be zero.
        time_edges = lr.get_n(int, rmfit=True)
        if time_edges:
            if f.data_type == 'ctime' or f.data_type == 'cspec':
                data_file.set_time_binning('By Edge Index', 'binned',
                                           np.array(time_edges))
        elif f.data_type == 'tte':
            # Read TI file
            if ti_file:
                with open(ti_file, 'r') as fp:
                    txt = list(fp)
                txt = txt[1:]
                tte_edges = np.array([t.strip() for t in txt], dtype=float)
                data_file.set_time_binning('By Time Edge', 'unbinned',
                                           np.array(tte_edges))
            else:
                warnings.warn("No TTE edges found.  Need '.ti' file")

        # time selections, if None is returned we need to read the next value anyway which should be zero.
        data_file.selections.source = transform_selections(
            lr.get_n(float, rmfit=True), drop_first=True)

        # background selections, if None is returned we need to read the next value anyway which should be zero.
        data_file.selections.background = transform_selections(
            lr.get_n(float, rmfit=True))

        # TODO: For now skip over binning names
        lr.skip(3)  # Assuming 'STACKED SPECTRA', 'LOG', 'LOG'

        # time and energy window ranges: (xmin, xmax, ymin, ymax)
        v = lr.get(4, float)
        data_file.views.time = View(v[0], v[1], v[2], v[3])
        v = lr.get(4, float)
        data_file.views.energy = View(v[0], v[1], v[2], v[3])

        # polynomial background order
        # data_file.background = {'poly_order': lr.get(cls=int)}
        poly_order = lr.get(cls=int)
        data_file.set_background_model('Polynomial', 'binned', poly_order)

        # Add the data file to a newly created lookup file
        lu = cls()
        lu.datafiles[data_file.filename] = data_file
        return lu

    def merge_lookup(self, lookup, overwrite=False):
        """Merge an existing lookup into this lookup

        Parameters:
        --------------
        lookup: GspecLookup
            The lookup object to be merged into this lookup
        overwrite: bool, optional
            If set to True, then any datanames in the current lookup will be overwritten
            if those same datanames are in the input lookup.  Default is False
        """
        # get datanames of the input lookup
        datanames = lookup.datafiles.keys()
        for dataname in datanames:
            # if dataname is already in this lookup and we don't want to overwrite
            if (dataname in self.datafiles) & (not overwrite):
                continue
            self.datafiles[dataname] = lookup.datafiles[dataname]

    # TODO: Remove?
    def split_off_dataname(self, dataname):
        """Return a new lookup object containing only the requested data file

        Parameters:
        --------------
        dataname: str
            The requested data filename

        Returns:
        -----------
        new_lookup: GspecLookup
            The new lookup object
        """
        self.assert_has_datafile(dataname)
        new_lookup = LookupFile()
        new_lookup[dataname] = self.datafiles[dataname]
        return new_lookup

    def display_lookup(self):
        """Pretty print a lookup for display (in json format)
        """
        lu = json.dumps(self.datafiles, indent=4, separators=(',', ': '),
                        cls=LookupEncoder)
        return lu


class LookupEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy arrays. Converts them to a list.
    """

    def default(self, obj):
        if isinstance(obj, DataFileLookup):
            d = dict(obj.__dict__)
            del d['filename']
            return d
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, Detector):
            return obj.short_name
        elif hasattr(obj, 'to_dict'):
            return obj.to_dict()
        elif hasattr(obj, '__dict__'):
            return obj.__dict__
        return json.JSONEncoder.default(self, obj)


class LookupDecoder(json.JSONDecoder):
    """Custom JSON decoder to turn JSON lists into numpy arrays
    """

    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook,
                                  *args, **kwargs)

    def object_hook(self, obj):
        # if object is a dictionary
        if type(obj) == dict:
            for key in obj.keys():
                # and if the value is a list, change to numpy array
                obj_type = type(obj[key])
                if obj_type == list:
                    obj[key] = np.array(obj[key], dtype=type(obj[key]))

        return obj
