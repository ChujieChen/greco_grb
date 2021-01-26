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

from unittest import TestCase
import datetime

from gbm.time import *
import csv


class TestTime(TestCase):
    def test_bn(self):
        # Test against historic BN numbers
        with open("bn_test.csv", "r") as bn_file:
            bn_reader = csv.DictReader(bn_file)
            for row in bn_reader:
                m = Met(float(row['trigger_time']))
                msg = "For met = %f\n" \
                      "Expected : '%s'\n" \
                      "Actual   : '%s'" % (m.met, row['burst_number'], m.bn)
                self.assertEqual(m.bn, row['burst_number'], msg)

    def test_unix_to_met(self):
        unix = 1329307200       # Wed, 15 Feb 2012 12:00:00 GMT
        expect = 351000002.000  # Wed, 15 Feb 2012 12:00:00 GMT

        m = Met.from_unix(unix)
        self.assertAlmostEqual(expect, m.met, places=6)

        unix = 1244628900       # Wed, 10 Jun 2009 10:15:00 GMT
        expect = 266321702.000  # Wed, 10 Jun 2009 10:15:00 GMT

        m = Met.from_unix(unix)
        self.assertAlmostEqual(expect, m.met, places=6)

        unix = 1221084000       # Wed, 10 Sep 2008 22:00:00 GMT
        expect = 242776801.000  # Wed, 10 Sep 2008 22:00:00 GMT

        m = Met.from_unix(unix)
        self.assertAlmostEqual(expect, m.met, places=6)

        unix = 1461164288       # Wed, 20 Apr 2016 14:58:08 GMT
        expect = 482857092.000  # Wed, 20 Apr 2016 14:58:08 GMT

        m = Met.from_unix(unix)
        self.assertAlmostEqual(expect, m.met, places=6)

    def test_utc_to_met(self):

        # Wed, 15 Feb 2012 12:00:00 GMT
        utc = datetime.datetime(2012, 2, 15, 12, 0)
        expect = 351000002.000

        m = Met.from_datetime(utc)
        self.assertAlmostEqual(expect, m.met, places=6)

        # Wed, 10 Jun 2009 10:15:00 GMT
        utc = datetime.datetime(2009, 6, 10, 10, 15)
        expect = 266321702.000

        m = Met.from_datetime(utc)
        self.assertAlmostEqual(expect, m.met, places=6)

        # Wed, 10 Sep 2008 22:00:00 GMT
        utc = datetime.datetime(2008, 9, 10, 22, 00)
        expect = 242776801.000

        m = Met.from_datetime(utc)
        self.assertAlmostEqual(expect, m.met, places=6)

        # Wed, 20 Apr 2016 14:58:08 GMT
        utc = datetime.datetime(2016, 4, 20, 14, 58, 8)
        expect = 482857092.000

        m = Met.from_datetime(utc)
        self.assertAlmostEqual(expect, m.met, places=6)

    def test_now(self):

        met = Met.now()
        print(met)

    def test_frac_of_day(self):
        d = datetime.datetime(2016, 1, 1, 12, 30)
        self.assertEqual(hms_to_fraction_of_day(d), 521)

    def check_leap_second(self, beg_next_day_utc, beg_next_day_met):
        # Test utc->met
        self.assertAlmostEqual(Met.from_unix(beg_next_day_utc - 2.0).met, beg_next_day_met - 3.0, places=6,
                               msg="Unix->Met m:58 seconds")
        self.assertAlmostEqual(Met.from_unix(beg_next_day_utc - 1.0).met, beg_next_day_met - 2.0, places=6,
                               msg="Unix->Met m:59 seconds")
        # (beg_next_day_met - 1) = No Unix 60 seconds
        self.assertAlmostEqual(Met.from_unix(beg_next_day_utc).met, beg_next_day_met, places=6,
                               msg="Unix->Met m+1:00 seconds")
        self.assertAlmostEqual(Met.from_unix(beg_next_day_utc + 1.0).met, beg_next_day_met + 1.0, places=6,
                               msg="Unix->Met m+1:01 seconds")

        # Test met->utc
        self.assertAlmostEqual(Met(beg_next_day_met - 3.0).unix, beg_next_day_utc - 2.0, places=6,
                               msg="Met->Unix m:58 seconds")
        self.assertAlmostEqual(Met(beg_next_day_met - 2.0).unix, beg_next_day_utc - 1.0, places=6,
                               msg="Met->Unix m:59 seconds")
        self.assertAlmostEqual(Met(beg_next_day_met - 1.0).unix, beg_next_day_utc - 1.0, places=6,
                               msg="Met->Unix m:60 seconds (repeat :59)")
        self.assertAlmostEqual(Met(beg_next_day_met).unix, beg_next_day_utc, places=6,
                               msg="Met->Unix m+1:00 seconds")
        self.assertAlmostEqual(Met(beg_next_day_met + 1.0).unix, beg_next_day_utc + 1.0, places=6,
                               msg="Met->Unix m+1:01 seconds")

    def test_2017_leap_second(self):
        self.check_leap_second(1483228800.0, 504921605.0)

    def test_2015_leap_second(self):
        self.check_leap_second(1435708800.0, 457401604.0)

    def test_2012_leap_second(self):
        self.check_leap_second(1341100800.0, 362793603.0)

    def test_2008_leap_second(self):
        self.check_leap_second(1230768000.0, 252460802.0)

    def test_2005_leap_second(self):
        self.check_leap_second(1136073600.0, 157766401.0)

    def test_met_to_gps(self):
        # Time values from HEASARC's xTime website were used with the exception of GPS time.
        # LIGO's converter webpage was used to convert UTC (from MET) to GPS time

        met = Met(157766410)  # 2006-01-01 00:00:09.000 UTC
        self.assertAlmostEqual(met.gps, 820108823, places=6)

        met = Met(252460802)  # 2009-01-01 00:00:00.000 UTC
        self.assertAlmostEqual(met.gps, 914803215, places=6)

        met = Met(362793605)  # 2012-07-01 00:00:02.000 UTC
        self.assertAlmostEqual(met.gps, 1025136018, places=6)

        met = Met(457401613)  # 2015-07-01 00:00:09.000 UTC
        self.assertAlmostEqual(met.gps, 1119744026, places=6)

        met = Met(506174405)  # 2017-01-15 12:00:00.000 UTC
        self.assertAlmostEqual(met.gps, 1168516818, places=6)

    def test_gps_to_met(self):
        # LIGO's converter webpage was used to convert UTC to GPS time
        # HEASARC's converter webpage was used to convert UTC to MET time

        met = Met.from_gps(825598814.0)  # 2006-03-05 13:00:00.000 UTC
        self.assertAlmostEqual(met.met, 163256401.0, places=6)

        met = Met.from_gps(934369815.0)  # 2009-08-15 11:10:00.000 UTC
        self.assertAlmostEqual(met.met, 272027402.0, places=6)

        met = Met.from_gps(1030610731.0)  # 2012-09-02 08:45:15.000 UTC
        self.assertAlmostEqual(met.met, 368268318.0, places=6)

        met = Met.from_gps(1113102199.0)  # 2015-04-15 03:03:03.000 UTC
        self.assertAlmostEqual(met.met, 450759786.0, places=6)

        met = Met.from_gps(1180827063.0)  # 2017-06-06 23:30:45.000 UTC
        self.assertAlmostEqual(met.met, 518484650.0, places=6)

        with self.assertRaises(Exception):
            Met.from_gps(662342412.0)

    def test_datetime_range_from(self):
        dates = hours_range_from(5, dt=datetime.datetime(2016, 8, 2, 2, 15))
        expect = [
            datetime.datetime(2016, 8, 1, 22, 0),
            datetime.datetime(2016, 8, 1, 23, 0),
            datetime.datetime(2016, 8, 2, 0, 0),
            datetime.datetime(2016, 8, 2, 1, 0),
            datetime.datetime(2016, 8, 2, 2, 0),
        ]
        count = 0
        for d in dates:
            self.assertEqual(d, expect[count])
            count += 1

    def test_date_range_from(self):
        dates = dates_range_from(5, dt=datetime.date(2016, 8, 2))
        expect = [
            datetime.date(2016, 7, 29),
            datetime.date(2016, 7, 30),
            datetime.date(2016, 7, 31),
            datetime.date(2016, 8, 1),
            datetime.date(2016, 8, 2),
        ]
        count = 0
        for d in dates:
            self.assertEqual(d, expect[count])
            count += 1

    def test_inclusive_date_range_inc(self):
        expect = [
            datetime.date(2016, 7, 29),
            datetime.date(2016, 7, 30),
            datetime.date(2016, 7, 31),
            datetime.date(2016, 8, 1),
            datetime.date(2016, 8, 2),
        ]
        dates = inclusive_date_range(expect[0], expect[-1])
        self.assertEqual(dates, expect)

    def test_inclusive_date_range_dec(self):
        expect = [
            datetime.date(2016, 8, 2),
            datetime.date(2016, 8, 1),
            datetime.date(2016, 7, 31),
            datetime.date(2016, 7, 30),
            datetime.date(2016, 7, 29),
        ]
        dates = inclusive_date_range(expect[0], expect[-1], datetime.timedelta(days=-1))
        self.assertEqual(dates, expect)

    def test_dates_from_hours(self):

        expect = [
            datetime.date(2016, 7, 29),
            datetime.date(2016, 7, 30),
            datetime.date(2016, 7, 31),
            datetime.date(2016, 8, 1),
            datetime.date(2016, 8, 2),
        ]

        hours = inclusive_date_range(datetime.datetime(2016, 7, 29, 0, 0), datetime.datetime(2016, 8, 2, 23, 0),
                                     step=datetime.timedelta(hours=1))
        dates = dates_from_hours(hours)
        self.assertEqual(dates, expect)


