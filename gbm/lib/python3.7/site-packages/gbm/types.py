# types.py: Module containing specialized types
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


def Struct(**kwargs):
    """Creates a Mutable struct that uses slots to discourage monkey patching as well as being more
        memory efficient.
    """

    class Struct:
        __slots__ = list(kwargs.keys())

        def __repr__(self):
            items = ("{}={!r}".format(k, getattr(self, k)) for k in
                     self.__slots__)
            return "{}({})".format(type(self).__name__, ", ".join(items))

        def __len__(self):
            return len(self.__slots__)

        def todict(self):
            d = dict()
            for k in self.__slots__:
                d[k] = getattr(self, k, None)
            return d

    s = Struct()
    for key, value in kwargs.items():
        setattr(s, key, value)
    return s


class ListReader:
    """Simple class that allows a list to be parsed without having to keep up with current position"""

    def __init__(self, in_list):
        self.pos = 0
        self.data = in_list

    def get(self, num_items=1, cls=None):
        result = None
        if num_items == 1:
            if cls:
                if cls == list:
                    result = [self.data[self.pos]]
                else:
                    result = cls(self.data[self.pos])
            else:
                result = self.data[self.pos]
        elif num_items > 1:
            if cls:
                result = list()
                for i in self.data[self.pos:self.pos + num_items]:
                    if cls == list:
                        result.append([i])
                    else:
                        result.append(cls(i))
            else:
                result = list(self.data[self.pos:self.pos + num_items])
        elif num_items < 0:
            raise ValueError("Not a valid number of requested items")
        self.pos += num_items
        return result

    def get_n(self, cls=None, rmfit=False):
        """Get a variable number of items based on the value at the current position"""
        num = self.get(cls=int)
        # Special case #1: RMFit had a bug where zero length still wrote a single zero value array.
        if num == 0 and rmfit:
            x = self.get(cls=int)
            if x != 0:
                raise ValueError(
                    "Non-zero value returned for array element in RMFit mode when zero was expected.")

        # Special case #2: If the array has only one element, it still needs to be returned in a list.
        if num == 1:
            return [self.get(num, cls)]

        return self.get(num, cls)

    def reset(self):
        self.pos = 0

    def seek(self, pos):
        self.pos = pos

    def skip(self, num):
        self.pos += num
        if self.pos > len(self.data):
            raise IndexError("Skipped past end of data")

    def rewind(self, num):
        self.pos -= num
        if self.pos < 0:
            raise IndexError("Rewound before beginning of data")
