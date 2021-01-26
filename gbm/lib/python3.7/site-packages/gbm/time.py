# time.py: Module containing time-related functions
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

import datetime
import warnings

from astropy.time import Time
from astropy.time.formats import TimeFromEpoch


class TimeFermiSec(TimeFromEpoch):
    """Represents the number of seconds elapsed since Jan 1, 2001 00:00:00 UTC 
    including leap seconds"""

    name = 'fermi'
    unit = 1.0 / 86400  # in days (1 day == 86400 seconds)
    epoch_val = '2001-01-01 00:01:04.184'
    epoch_val2 = None
    epoch_scale = 'tt'  # Scale for epoch_val class attribute
    epoch_format = 'iso'  # Format for epoch_val class attribute


def round_half_to_nearest_even(num):
    n = int(num)
    v = abs(num - n)
    if v > 0.5 or (v == 0.5 and n % 2):
        return n + 1 - 2 * int(n < 0)
    else:
        return n


def hms_to_fraction_of_day(value):
    """The fraction of day as computed by the original ops code.

    Args:
        value (:class:`datetime.datetime`): The date/time
    
    Returns:
        float: The fraction of day
    """

    result = round_half_to_nearest_even(((
                                                     value.hour * 3600 + value.minute * 60 + value.second) / 86400) * 1000)
    return min(result, 999)


def fraction_of_day_to_hms(value):
    """The hour, minute, second for a given fraction of day
    
    Args:
        value (float): The fraction of day
        
    Returns:
        (int, int, int): Hour, minute, second
    """
    s = int((value / 1000) * 86400)
    h = s // 3600
    s -= h * 3600
    m = s // 60
    s -= m * 60
    return h, m, s


class Met:
    """Class representing the Fermi MET Epoch and allowing time conversions
    to and from it.
    
    Parameters:
        secs (float): The MET
    
    Attributes:
        bn (str): The MET converted to bust number format: 'YYMMDDfff'
        datetime (:class:`datetime.datetime`): A datetime object for the MET
        gps (float): The number of seconds since Jan 6, 1980 00:00:00 
                    (leap seconds are removed)
        jd (float): The Julian Date associated with the MET
        met (float): The MET
        mjd (float): The modified Julian Date associated with the MET
        time (:class:`astropy.time.Time`): The astropy time object for the MET
        unix (float): The number of seconds since Jan 1, 1970 00:00:00 with 
                      the leap seconds removed
        ymd (str): The MET converted to the form `YYMMDD` in UTC
        ymd_h (str): The MET converted to the form of YYMMDD_HHz in UTC
    """

    # Mission Elapsed Time (Number of seconds since 2001-01-01 00:00:00 UTC)

    def __init__(self, secs):
        """Creates a Met object with the time set to the number of seconds since Jan 1, 2001 00:00:00 UTC including the
         leap seconds"""
        if secs < 0:
            warnings.warn("Time before GBM mission epoch")
            # raise Exception("Time before GBM mission epoch")
        self.__time = Time(secs, format='fermi')

    @classmethod
    def from_iso(cls, str_time):
        """Create a new Met object from an ISO-format UTC string
    
        Args:
            str_time (str): The ISO string
        
        Returns:
            :class:`Met`: The Met object
        """
        if '.' in str_time:
            dt = datetime.datetime.strptime(str_time, '%Y-%m-%dT%H:%M:%S.%f')
        else:
            dt = datetime.datetime.strptime(str_time, '%Y-%m-%dT%H:%M:%S')
        return cls.from_datetime(dt)

    @property
    def met(self):
        return self.__time.fermi

    # Astropy Time
    @property
    def time(self):
        return self.__time

    @classmethod
    def from_time(cls, atime):
        """Creates a new Met object from an astropy.Time object
    
        Args:
            atime (:class:`astropy.time.Time`): The astropy time object
        
        Returns:
            :class:`Met`: The Met object
        """
        obj = cls(0)
        obj.__time = atime
        if obj.met < 0:
            raise Exception("Time before GBM mission epoch")
        return obj

    # Python's datetime
    @property
    def datetime(self):
        try:
            return self.__time.utc.to_datetime(datetime.timezone.utc)
        except ValueError:
            # Repeat last met for a leap second
            return Met(self.met - 1).datetime

    @classmethod
    def from_datetime(cls, dt):
        """Creates a new Met object from a datetime.datetime object
    
        Args:
            dt (:class:`datetime.datetime`): The datetime object
        
        Returns:
            :class:`Met`: The Met object
        """
        return cls.from_time(Time(dt, format='datetime'))

    # Unix timestamp (Number of seconds since 1970-01-01 00:00:00 UTC (leap seconds are ignored))
    @property
    def unix(self):
        return self.datetime.timestamp()

    @classmethod
    def from_unix(cls, unix):
        """Creates a new Met object from a Unix timestamp
    
        Args:
            unix (float): A Unix time
        
        Returns:
            :class:`Met`: The Met object
        """
        return cls.from_datetime(datetime.datetime.utcfromtimestamp(unix))

    # GPS timestamp (Number of seconds since Jan 6, 1980 00:00:00 UTC (leap seconds are ignored))
    @property
    def gps(self):
        return self.__time.gps

    @classmethod
    def from_gps(cls, gps):
        """Creates a new Met object from a GPS timestamp
    
        Args:
            gsp (float): A GPS time
        
        Returns:
            :class:`Met`: The Met object
        """
        return cls.from_time(Time(gps, format='gps'))

    # Julian date
    @property
    def jd(self):
        return self.__time.jd

    @classmethod
    def from_jd(cls, jd):
        """Creates a new Met object from a Julian Date
    
        Args:
            jd (float): A Julian Date
        
        Returns:
            :class:`Met`: The Met object
        """
        return cls.from_time(Time(jd, format='jd'))

    # Modified Julian Date
    @property
    def mjd(self):
        return self.__time.utc.mjd

    @classmethod
    def from_mjd(cls, mjd):
        """Creates a new Met object from a Modified Julian Date
    
        Args:
            mjd (float): A Modified Julian Date
        
        Returns:
            :class:`Met`: The Met object
        """
        return cls.from_time(Time(mjd, format='mjd'))

    # GBM Burst Number (YYMMDDFFF)
    @property
    def bn(self):

        # Adjust to match a known bug in the old pipeline software
        adj_met = self.met
        if 157766399.0 < adj_met < 252460801.0:
            adj_met += 1
        elif 252460800.0 < adj_met <= 253497600.0:
            adj_met += 2

        # To ensure compatibility with the number produced by the pipeline, we are doing it the inefficient way
        m = Met(adj_met)
        utc_val = m.datetime
        fraction = hms_to_fraction_of_day(utc_val)

        return "{}{:03d}".format(utc_val.strftime("%y%m%d"), fraction)

    @classmethod
    def from_bn(cls, bn):
        """Creates a new Met object from a 'YYMMDDfff' string
    
        Args:
            bn (str): A burst number string
        
        Returns:
            :class:`Met`: The Met object
        """
        dt = datetime.datetime.strptime(bn[:6], '%y%m%d')
        hms = fraction_of_day_to_hms(int(bn[6:]))
        dt = datetime.datetime(dt.year, dt.month, dt.day, hms[0], hms[1],
                               hms[2], tzinfo=datetime.timezone.utc)
        obj = cls.from_datetime(dt)

        # Adjust to match a known bug in the old pipeline software
        adj_met = obj.met
        if 157766400.0 < adj_met < 252460802.0:
            adj_met -= 1
        elif 252460802.0 < adj_met <= 253497602.0:
            adj_met -= 2

        return Met(adj_met)

    # Year, Month, and Day as YYMMDD
    @property
    def ymd(self):
        return self.datetime.strftime("%y%m%d")

    @classmethod
    def from_ymd(cls, ymd):
        """Creates a new Met object from a 'YYMMDD' string
    
        Args:
            ymd (str): A YYMMDD string
        
        Returns:
            :class:`Met`: The Met object
        """
        dt = datetime.datetime.strptime(ymd, '%y%m%d')
        return cls.from_datetime(dt)

    # Year, Month, Day, and Hour as YYMMDD_HH
    @property
    def ymd_h(self):
        return self.datetime.strftime("%y%m%d_%Hz")

    @classmethod
    def from_ymd_h(cls, ymd):
        """Creates a new Met object from a 'YYMMDD_HHz' string
    
        Args:
            ymd (str): A YYMMDD_HHz string
        
        Returns:
            :class:`Met`: The Met object
        """
        dt = datetime.datetime.strptime(ymd, '%y%m%d_%Hz')
        return cls.from_datetime(dt)

    # Current time
    @classmethod
    def now(cls):
        """Creates a new Met object from the current time
    
        Returns:
            :class:`Met`: The Met object
        """
        m = cls(0)
        m.__time = Time.now()
        return m

    # String functions

    def iso(self):
        """Returns the MET value as a string in the form of 
        yyyy-mm-ddTHH:MM:SS in UT
    
        Returns:
            :str: the ISO string
        """
        return self.datetime.strftime("%Y-%m-%dT%H:%M:%S")

    def __repr__(self):
        """Returns a string representation of the Met object"""
        return "<Met seconds = {:.6f}>".format(self.met)

    # Math functions
    def add(self, x):
        """Returns an Met object with its value set to this object's value 
        with x seconds added to it. Can also use the ``+`` operator.
    
        Args:
            x (float): seconds to add
        
        Returns:
            :class:`Met`: The Met object
        """
        if not (isinstance(x, int) or isinstance(x, float)):
            raise ValueError("Can only add int or float to Met")
        return Met(self.met + x)

    def sub(self, x):
        """Returns an Met object with its value set to this object's value 
        with x seconds subtracted from it. Can also use the ``-`` operator.
    
        Args:
            x (float): seconds to subtract
        
        Returns:
            :class:`Met`: The Met object
        """
        if isinstance(x, Met):
            return self.met - x.met
        elif isinstance(x, int) or isinstance(x, float):
            return Met(self.met - x)
        raise ValueError("Can only subtract int, float or Met from Met")

    # Overriding built-in operators
    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.sub(other)

    def __lt__(self, other):
        if isinstance(other, Met):
            return self.met < other.met
        else:
            raise TypeError(
                "'<' not supported between instances of 'Met' and '{}'".format(
                    type(other)))

    def __le__(self, other):
        if isinstance(other, Met):
            return self.met <= other.met
        else:
            raise TypeError(
                "'<=' not supported between instances of 'Met' and '{}'".format(
                    type(other)))

    def __gt__(self, other):
        if isinstance(other, Met):
            return self.met > other.met
        else:
            raise TypeError(
                "'>' not supported between instances of 'Met' and '{}'".format(
                    type(other)))

    def __ge__(self, other):
        if isinstance(other, Met):
            return self.met >= other.met
        else:
            raise TypeError(
                "'>=' not supported between instances of 'Met' and '{}'".format(
                    type(other)))

    def __eq__(self, other):
        if isinstance(other, Met):
            return self.met == other.met
        else:
            raise TypeError(
                "'==' not supported between instances of 'Met' and '{}'".format(
                    type(other)))

    def __ne__(self, other):
        if isinstance(other, Met):
            return self.met != other.met
        else:
            raise TypeError(
                "'!=' not supported between instances of 'Met' and '{}'".format(
                    type(other)))


# Some time related functions

def inclusive_date_range(start, stop, step=datetime.timedelta(days=1)):
    """Creates a list of Met from start to stop times
    
    Args:
        start (:class:`Met`): The start MET
        stop (:class:`Met`): The end MET
        step (:class:`datetime.timedelta, optional): 
            The step size. Default is 1 day.
    
    Returns:
        list of :class:`Met`: The list of Met objects
    """
    d = start
    result = []

    if start <= stop:
        earliest, latest = start, stop
    else:
        earliest, latest = stop, start

    while earliest <= d <= latest:
        result.append(d)
        d += step

    return result


def dates_range_from(num_days, dt=datetime.datetime.utcnow().date()):
    """Creates a list of dates within the given range
    
    Args:
        num_days (int): Number of days to include in the list
        dt (:class:`datetime.date`, optional): 
            The last date to be included in the list. Default is current date.
    
    Returns:
        list: List of date values representing hours.
    """
    d = dt - datetime.timedelta(days=num_days - 1)
    return inclusive_date_range(d, dt)


def hours_range_from(num_hours, dt=datetime.datetime.utcnow()):
    """Creates a list of datetimes within the given range
    
    Args:
        num_hours (int): Number of hours to include in the list
        dt (:class:`datetime.datetime`, optional): 
            The last hour to be included in the list (datetime will be 
            truncated to hour value). Default is current hour.
    
    Returns:
        list: List of datetime values representing hours.
    """
    d = datetime.datetime(dt.year, dt.month, dt.day, dt.hour, 0, 0)
    d -= datetime.timedelta(hours=num_hours - 1)

    return inclusive_date_range(d, dt, datetime.timedelta(hours=1))


def dates_from_hours(hours):
    """Converts a list of hours to a list of days spanned
    
    Args:
        hours (list of :class:`datetime.date`): List of hours
    
    Returns:
        list: The list of dates
    """
    return inclusive_date_range(hours[0].date(), hours[-1].date())
