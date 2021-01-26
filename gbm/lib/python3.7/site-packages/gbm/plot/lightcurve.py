# lightcurve.py: Plot class for lightcurves
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

from .gbmplot import GbmPlot, Histo, HistoErrorbars, HistoFilled, \
    LightcurveBackground
from .lib import *


class Lightcurve(GbmPlot):
    """Class for plotting lightcurves and lightcurve paraphernalia.
    
    Parameters:
        data (:class:`~gbm.data.primitives.TimeBins`, optional): 
            The lightcurve data to plot
        background (:class:`~gbm.background.BackgroundRates`, optional): 
            The background rates to plot
        **kwargs: Options to pass to :class:`~.gbmplot.GbmPlot`
    
    Attributes:
        ax (:class:`matplotlib.axes`): The matplotlib axes object for the plot
        background (:class:`~.gbmplot.LightcurveBackground`): 
            The background plot element
        canvas (Canvas Backend object): The plotting canvas, if set upon 
                                        initialization.
        errorbars (:class:`~.gbmplot.HistoErrorbars`): The error bars plot element
        fig (:class:`matplotlib.figure`): The matplotlib figure object
        lightcurve (:class:`~.gbmplot.Histo`): The lightcurve plot element
        selections (list of :class:`~.gbmplot.HistoFilled`):
            The list of selection plot elements
        xlim (float, float): The plotting range of the x axis. 
                             This attribute can be set.
        xscale (str): The scale of the x axis, either 'linear' or 'log'. 
                      This attribute can be set.
        ylim (float, float): The plotting range of the y axis. 
                             This attribute can be set.
        yscale (str): The scale of the y axis, either 'linear' or 'log'. 
                      This attribute can be set.
    """

    def __init__(self, data=None, background=None, canvas=None, axis=None,
                 **kwargs):
        super().__init__(canvas=canvas, axis=axis, **kwargs)

        self._lc = None
        self._errorbars = None
        self._bkgd = None
        self._selections = []

        # initialize the plot axes, labels, ticks, and scales
        self._ax.set_xlabel('Time (s)', fontsize=PLOTFONTSIZE)
        self._ax.set_ylabel('Count Rate (count/s)', fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('linear')
        self._ax.set_yscale('linear')

        # plot data and/or background if set on init
        if data is not None:
            self.set_data(data)
            self._ax.set_xlim(data.range)

            self._ax.set_ylim(0.9 * np.min(data.rates),
                              1.1 * np.max(data.rates))
        if background is not None:
            self.set_background(background)

    @property
    def lightcurve(self):
        return self._lc

    @property
    def errorbars(self):
        return self._errorbars

    @property
    def background(self):
        return self._bkgd

    @property
    def selections(self):
        return self._selections

    def set_data(self, data):
        """Set the lightcurve plotting data. If a lightcurve already exists,
        this triggers a replot of the lightcurve.
        
        Args:
            data (:class:`~gbm.data.primitives.TimeBins`): 
                The lightcurve data to plot
        """
        lc_color, lc_alpha, lc_kwargs = self._lc_settings()
        self._lc = Histo(data, self._ax, color=lc_color, alpha=lc_alpha,
                         **lc_kwargs)
        eb_color, eb_alpha, eb_kwargs = self._eb_settings()
        self._errorbars = HistoErrorbars(data, self._ax, color=eb_color,
                                         alpha=eb_alpha, **eb_kwargs)

    def add_selection(self, data):
        """Add a selection to the plot.  This adds a new selection to a list
        of existing selections.
        
        Args:
            data (:class:`~gbm.data.primitives.TimeBins`): 
                The lightcurve data selection to plot
        """
        color, alpha, kwargs = self._selection_settings()
        select = HistoFilled(data, self._ax, color=color, alpha=alpha,
                             **kwargs)
        self._selections.append(select)

    def set_background(self, background):
        """Set the background plotting data. If a background already exists,
        this triggers a replot of the background.
        
        Args:
            background (:class:`~gbm.background.BackgroundRates`): 
                The background model to plot
        """
        color, cent_color, err_color, alpha, cent_alpha, err_alpha, \
            kwargs = self._bkgd_settings()
        self._bkgd = LightcurveBackground(background, self._ax, color=color,
                                          cent_color=cent_color,
                                          err_color=err_color,
                                          alpha=alpha, cent_alpha=BKGD_ALPHA,
                                          err_alpha=BKGD_ERROR_ALPHA,
                                          zorder=1000)

    def remove_data(self):
        """Remove the lightcurve from the plot.
        """
        self._lc.remove()
        self._lc = None

    def remove_errorbars(self):
        """Remove the lightcurve error bars from the plot.
        """
        self._errorbars.remove()
        self._errorbars = None

    def remove_background(self):
        """Remove the background from the plot.
        """
        self._bkgd.remove()
        self._bkgd = None

    def remove_selections(self):
        """Remove the selections from the plot.
        """
        [selection.remove() for selection in self._selections]
        self._selections = []

    def _lc_settings(self):
        """The default settings for the lightcurve. If a lightcurve already
        exists, use its settings instead.
        """
        if self._lc is None:
            lc_color = DATA_COLOR
            lc_alpha = None
            lc_kwargs = {}
        else:
            lc_color = self._lc.color
            lc_alpha = self._lc.alpha
            lc_kwargs = self._lc._kwargs
        return (lc_color, lc_alpha, lc_kwargs)

    def _eb_settings(self):
        """The default settings for the errorbars. If a lightcurve already
        exists, use its errorbars settings instead.
        """
        if self._errorbars is None:
            eb_color = DATA_ERROR_COLOR
            eb_alpha = None
            eb_kwargs = {}
        else:
            eb_color = self._errorbars.color
            eb_alpha = self._errorbars.alpha
            eb_kwargs = self._errorbars._kwargs
        return (eb_color, eb_alpha, eb_kwargs)

    def _bkgd_settings(self):
        """The default settings for the background. If a background already
        exists, use its settings instead.
        """
        if self._bkgd is None:
            color = BKGD_COLOR
            cent_color = None
            err_color = None
            alpha = None
            cent_alpha = BKGD_ALPHA
            err_alpha = BKGD_ERROR_ALPHA
            kwargs = {'linewidth': BKGD_WIDTH}
        else:
            color = self._bkgd.color
            cent_color = self._bkgd.cent_color
            err_color = self._bkgd.err_color
            alpha = self._bkgd.alpha
            cent_alpha = self._bkgd.cent_alpha
            err_alpha = self._bkgd.err_alpha
            kwargs = self._bkgd._kwargs
        return color, cent_color, err_color, alpha, cent_alpha, err_alpha, kwargs

    def _selection_settings(self):
        """The default settings for a selection. If a selection already
        exists, use its settings instead.
        """
        if len(self._selections) == 0:
            color = DATA_SELECTED_COLOR
            alpha = DATA_SELECTED_ALPHA
            kwargs = {}
        else:
            color = self._selections[0].color
            alpha = self._selections[0].alpha
            kwargs = self._selections[0]._kwargs
        return color, alpha, kwargs
