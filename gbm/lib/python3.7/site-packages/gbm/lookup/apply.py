# apply.py: Module contain functions to apply actions from lookups
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
from ..data.phaii import PHAII, TTE
from ..background.background import BackgroundFitter


def rebin_time(data_obj, data_lookup):
    """Temporal rebin based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    phaii: PHAII
        The rebinned PHAII object
    """
    if isinstance(data_obj, TTE):
        func = data_obj.to_phaii
    elif isinstance(data_obj, PHAII):
        func = data_obj.rebin_time
    else:
        raise TypeError('Data must either be PHAII or TTE')

    binnings = data_lookup.binnings.time
    if binnings is None:
        return data_obj

    for binning in binnings:
        phaii = func(binning.method, *binning.args, time_range=(binning.start,
                                                                binning.stop))
        func = phaii.rebin_time
    return phaii


def rebin_energy(data_obj, data_lookup):
    """Energy rebin based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    phaii: PHAII
        The rebinned PHAII object
    """

    try:
        func = data_obj.rebin_energy
    except:
        raise TypeError('Not a valid data object for energy rebinning')

    binnings = data_lookup.binnings.energy
    if binnings is None:
        return data_obj

    for binning in binnings:
        phaii = func(binning.method, *binning.args, emin=binning.start,
                     emax=binning.stop)
        func = phaii.rebin_time
    return phaii


def source_selection(data_obj, data_lookup):
    """Source selection (temporal) based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    new_obj: PHAII or TTE
        A new data object containing the source selection
    """
    try:
        func = data_obj.slice_time
    except:
        raise TypeError('Not a valid data object for source selection')

    new_obj = func(data_lookup.selections.source)
    return new_obj


def energy_selection(data_obj, data_lookup):
    """Energy selection based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    new_obj: PHAII or TTE
        A new data object containing the energy selection
    """
    try:
        func = data_obj.slice_energy
    except:
        raise TypeError('Not a valid data object for energy selection')

    new_obj = func(data_lookup.selections.energy)
    return new_obj


def background_selection(data_obj, data_lookup):
    """Background selection (temporal) based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    new_obj: PHAII or TTE
        A new data object containing the background selection
    """
    try:
        func = data_obj.slice_time
    except:
        raise TypeError('Not a valid data object for background selection')

    new_obj = func(data_lookup.selections.background)
    return new_obj


def background_fit(data_obj, data_lookup):
    """Perform a fit of the background based on a lookup
    
    Parameters:
    -----------
    data_obj: PHAII or TTE
        The data to rebin
    data_lookup: DataFileLookup
        The lookup for the data
    
    Returns:
    -----------
    fitter: BackgroundFitter
        The fitter class containing the fitted background model
    """
    if isinstance(data_obj, PHAII):
        fitter_func = BackgroundFitter.from_phaii
    elif isinstance(data_obj, TTE):
        fitter_func = BackgroundFitter.from_tte
    else:
        raise TypeError('Not a valid data object for background fitting')

    fitter = fitter_func(data_obj, data_lookup.background.method,
                         time_ranges=data_lookup.selections.background)
    fitter.fit(*data_lookup.background.args, **data_lookup.background.kwargs)

    return fitter


def lightcurve_view(plot_obj, data_lookup):
    """Set the view range of the lightcurve based on a lookup
    
    Parameters:
    -----------
    plot_obj: Lightcurve
        The lightcurve plot object
    data_lookup: DataFileLookup
        The lookup for the data
    """
    try:
        plot_obj.xlim = data_lookup.views.time.xrange()
        plot_obj.ylim = data_lookup.views.time.yrange()
    except:
        pass


def spectrum_view(plot_obj, data_lookup):
    """Set the view range of the spectrum based on a lookup
    
    Parameters:
    -----------
    plot_obj: Spectrum
        The spectrum plot object
    data_lookup: DataFileLookup
        The lookup for the data
    """
    try:
        plot_obj.xlim = data_lookup.views.energy.xrange()
        plot_obj.ylim = data_lookup.views.energy.yrange()
    except:
        pass
