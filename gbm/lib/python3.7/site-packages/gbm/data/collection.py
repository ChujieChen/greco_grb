# collection.py: Data collection classes
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
import re
from functools import partial
from collections import OrderedDict


class DataCollection:
    """A container for a collection of like data objects, such as a collection
    of :class:`Ctime` objects.  This class exposes the individual objects'
    attributes, and it exposes the methods of the object class so that methods 
    can be called on the collection as a whole.  For that reason, each object in 
    the DataCollection must be of the same type, otherwise an error is raised.  
    The type of the collection is set by the first object inserted into the 
    collection and the collection type is immutable thereafter.
    
    Objects are stored in the collection in the order they are added.
    
    The number of items in the collection can be retrieved by ``len()`` and
    one can iterate over the items:: 
        [data_item for data_item in DataCollection]
    
    In addition to the DataCollection methods, all of the individual object 
    attributes and methods are exposed, and they become methods of the 
    DataCollection.  Note that individual object attributes become *methods* 
    i.e. if you have an item attribute called item.name, then the corresponding 
    DataCollection method would be item.name().
    
    Attributes:
        items (list): The names of the items in the DataCollection
        types (str): The type of the objects in the DataCollection            
    """

    def __init__(self):
        self._data_dict = OrderedDict()
        self._type = None

    def __iter__(self):
        for item in self._data_dict.values():
            yield item

    def __len__(self):
        return len(self._data_dict)

    def _enforce_type(self, data_item):
        if not isinstance(data_item, self._type) and self._type is not None:
            raise TypeError(
                'Incorrect data item for {}'.format(self.__class__.__name__))

    @property
    def items(self):
        return list(self._data_dict.keys())

    @property
    def types(self):
        return self._type

    @classmethod
    def from_list(cls, data_list, names=None):
        """Given a list of objects and optionally a list of corresponding names, 
        create a new DataCollection. 
        
        Args:
            data_list (list of :obj:`objects`): 
                The list of objects to be in the collection
            names (list of :obj:`str`, optional):  
                The list of corresponding names to the objects.  If not set, 
                will try to retrieve a name from object.filename (assuming it's 
                a data object). If that fails, each item will be named 
                ambiguously 'item1', 'item2', etc.
        
        Returns                
            :py:class:`DataCollection`: The newly created collection
        """
        obj = cls()

        # set the names
        if names is not None:
            if len(names) != len(data_list):
                raise ValueError('Names list must be same size as data list')
        else:
            names = [None] * len(data_list)

        # include the objects
        for data_item, name in zip(data_list, names):
            obj.include(data_item, name=name)

        return obj

    def to_list(self):
        """Return the objects contained in the DataCollection as a list.
        
        Returns:
            (list of :obj:`objects`): 
                The list of objects, in the order that they were inserted
        """
        return [self.get_item(name) for name in self.items]

    def include(self, data_item, name=None):
        """Insert an object into the collection.  The first item inserted will 
        set the immutable type.
        
        Args:
            data_item (:obj:`object`): A data object to include
            name (str, optional): 
                An optional corresponding name.  If not set, will try to 
                retrieve a name from object.filename (assuming it's a data 
                object). If that fails, each item will be named ambiguously 
                'item1', 'item2', etc.
        """
        # if this is the first item inserted, set the type of the Collection
        # and expose the attributes and methods of the object
        if len(self) == 0:
            self._type = type(data_item)
            dir = [key for key in data_item.__dir__() if
                   not re.match('_.', key)]
            for key in dir:
                setattr(self, key, partial(self._method_call, key))   
        else:
            # otherwise, ensure that each object inserted is of the same type
            if type(data_item) != self._type:
                raise TypeError('A DataCollection must contain like objects')

        # insert with user-defined name
        if name is not None:
            self._data_dict[name] = data_item
        else:
            # or try to insert using filename attribute
            try:
                self._data_dict[data_item.filename] = data_item
            # otherwise default to ambiguity
            except AttributeError:
                self._data_dict['item{}'.format(len(self) + 1)] = data_item

    def remove(self, item_name):
        """Remove an object from the collection given the name 
        
        Args:
            item_name (str): The name of the item to remove
        """
        self._data_dict.pop(item_name)

    def get_item(self, item_name):
        """Retrieve an object from the DataCollection by name
        
        Args:
            item_name (str): The name of the item to retrieve
        
        Returns:
            :obj:`object`: The retrieved data item
        """
        return self._data_dict[item_name]

    def _method_call(self, method_name, *args, **kwargs):
        """This is the wrapper for the exposde attribute and method calls.  
        Applies method_name over all items in the DataCollection
        
        Args:
            method_name (str): The name of the method or attribute
            *args: Additional arguments to be passed to the method
            **kwargs: Additional keyword arguments to be passed to the method
        
        Returns:
            None or list: If not None, will return the results from all 
            objects in the list     
       """
        # get the attributes/methods for each item
        refs = [getattr(obj, method_name) for obj in self._data_dict.values()]

        # if method_name is a method, then it will be callable
        if callable(refs[0]):
            res = [getattr(obj, method_name)(*args, **kwargs)
                   for obj in self._data_dict.values()]
        # otherwise, method_name will not be callable if it is an attribute
        else:
            # we are setting an attribute    
            if len(args) != 0:
                res = [setattr(obj, method_name, *args)
                       for obj in self._data_dict.values()]
            # we are retrieving an attribute
            else:
                res = refs

        if res[0] is not None:
            return res


class GbmDetectorCollection(DataCollection):
    """A container for a collection of GBM-specific data objects, such as a 
    collection of ``Ctime`` objects from different detectors.
    
    The special behavior of this class is to provide a way to interact with
    a collection of detector data that may contain a mix of different *types*
    of detectors.  For example, many times we want a collection of GBM NaI
    and GBM BGO detectors.  These detectors have very different energy ranges,
    and so may require different inputs for a variety of functions.  This
    collection allows one to specify the different arguments for NaI and BGO
    data without having to implement many ugly and space-wasting loops and
    ``if...else`` decisions.
    
    In addition to the GbmDetectorCollection methods, all of the individual 
    object attributes and methods are exposed, and they become methods of the 
    DataCollection.  Note that individual object attributes become *methods* 
    i.e. if you have an item attribute called item.name, then the corresponding 
    DataCollection method would be item.name().    

    Attributes:
        items (list): The names of the items in the DataCollection
        types (str): The type of the objects in the DataCollection            
    """

    def __init__(self):
        super().__init__()
        self._dets = []

    @classmethod
    def from_list(cls, data_list, names=None, dets=None):
        """Given a list of objects and optionally a list of corresponding names
        and corresponding detector names, create a new GbmDetectorCollection. 
        
        Args:
            data_list (list of :obj:`objects`): 
                The list of objects to be in the collection
            names (list of :obj:`str`, optional):  
                The list of corresponding names to the objects.  If not set, 
                will try to retrieve a name from object.filename (assuming it's 
                a data object). If that fails, each item will be named 
                ambiguously 'item1', 'item2', etc.
            dets (list of :obj:`str`, optional): 
                The detector names for each object. If not set, will try to
                retrieve from the object.detector attribute.  If that attribute 
                doesn't exist, an error will be raised, and the user will need 
                to specify this list.
        
        Returns                
            :py:class:`GbmDetectorCollection`: The newly created collection
        """
        obj = cls()

        # set the detector names
        if dets is not None:
            if len(dets) != len(data_list):
                raise ValueError(
                    'Detector list must be same size as data list')
        else:
            try:
                dets = [data_item.detector for data_item in data_list]
            except:
                raise AttributeError('Cannot find detector information. '
                                     'Need to manually set')
        # set the names
        if names is not None:
            if len(names) != len(data_list):
                raise ValueError('Names list must be same size as data list')
        else:
            names = [None] * len(data_list)

        # include the objects
        [obj.include(data_item, det, name=name) for (data_item, det, name)
         in zip(data_list, dets, names)]

        return obj

    def remove(self, item_name):
        """Remove an object from the collection given the name 
        
        Args:
            item_name (str): The name of the item to remove
        """
        index = [item == item_name for item in self.items].index(True)
        self._dets.pop(index)
        self._data_dict.pop(item_name)

    def include(self, data_item, det, name=None):
        """Insert an object into the GbmDetectorCollection.  The first item 
        inserted will set the immutable type.
        
        Args:
            data_item (:obj:`object`): A data object to include
            det (str): The corresponding detector for the item
            name (str, optional): 
                An optional corresponding name.  If not set, will try to 
                retrieve a name from object.filename (assuming it's a data 
                object). If that fails, each item will be named ambiguously 
                'item1', 'item2', etc.
        """
        super().include(data_item, name=None)
        self._dets.append(det)

    def _method_call(self, method_name, *args, nai_args=(), nai_kwargs=None,
                     bgo_args=(), bgo_kwargs=None, **kwargs):
        """This is the wrapper for the attribute and method calls.  Applies 
        method_name over all items in the GbmDetectorCollection.

        Args:
            method_name (str): The name of the method or attribute
            *args: Additional arguments to be passed to the method
            nai_args: Arguments to be applied only to the NaI objects
            bgo_args: Arguments to be applied only to the BGO objects
            nai_kwargs: Keywords to be applied only to the NaI objects
            bgo_kwargs: Keywords to be applied only to the BGO objects
            **kwargs: Additional keyword arguments to be passed to the 
                      method. Will be applied to both NaI and BGO objects 
                      and will be appended to any existing keywords from 
                      nai_kwargs or bgo_kwargs
        
        Returns:       
        None or list: If not None, will return the results from all objects 
                      in the list     
        """
        if nai_kwargs is None:
            nai_kwargs = {}
        if bgo_kwargs is None:
            bgo_kwargs = {}

        if len(args) > 0:
            nai_args = args
            bgo_args = args
        if len(kwargs) > 0:
            nai_kwargs.update(kwargs)
            bgo_kwargs.update(kwargs)

        # get the attributes/methods for each item
        refs = [getattr(obj, method_name) for obj in self._data_dict.values()]

        # if method_name is a method, then it will be callable
        if callable(refs[0]):
            res = []
            for obj, det in zip(self._data_dict.values(), self._dets):
                if 'n' in det:
                    our_args = nai_args
                    our_kwargs = nai_kwargs
                elif 'b' in det:
                    our_args = bgo_args
                    our_kwargs = bgo_kwargs
                res.append(getattr(obj, method_name)(*our_args, **our_kwargs))

        # otherwise, method_name will not be callable if it is an attribute
        else:
            # we are setting an attribute    
            if len(nai_args) != 0 or len(bgo_args) != 0:
                res = []
                for obj, det in zip(self._data_dict.values(), self._dets):
                    if 'n' in obj.detector:
                        our_args = nai_args
                    elif 'b' in obj.detector:
                        our_args = bgo_args
                    res.append(setattr(obj, method_name, *args))
            # we are retrieving an attribute
            else:
                res = refs

        if res[0] is not None:
            return res
