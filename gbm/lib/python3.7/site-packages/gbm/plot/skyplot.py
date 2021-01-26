# skyplot.py: Plot class for Fermi observing sky maps
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
from .gbmplot import Collection
from .gbmplot import DetectorPointing, GalacticPlane, SkyHeatmap, SkyPolygon
from .gbmplot import GbmPlot, SkyLine, SkyCircle, Sun
from .lib import *
from ..coords import get_sun_loc


class SkyPlot(GbmPlot):
    """Class for plotting on the sky in equatorial coordinates

    Parameters:
        projection (str, optional): The projection of the map. 
                                    Default is 'mollweide'
        flipped (bool, optional): 
            If True, the RA axis is flipped, following astronomical convention. \
            Default is True.
        **kwargs: Options to pass to :class:`~.gbmplot.GbmPlot`

    Attributes:
        ax (:class:`matplotlib.axes`): The matplotlib axes object for the plot
        background_color (str): The color of the plot background. This attribute
                                can be set.
        canvas (Canvas Backend object): The plotting canvas, if set upon 
                                        initialization.
        canvas_color (str): The color of the plotting canvas. This attribute 
                            can be set.
        detectors (:class:`~.gbmplot.Collection` of :class:`~.gbmplot.DetectorPointing`):
            The collection of detector plot elements
        earth (:class:`~.gbmplot.SkyCircle`): The Earth plot element
        fig (:class:`matplotlib.figure`): The matplotlib figure object
        fontsize (int): The font size of the text labels. This attribute can be set.
        gc (:class:`~.gbmplot.GalacticPlane`):
            The reference to the galactic plane plot element
        loc_contours: (:class:`~.gbmplot.Collection` of :class:`~.gbmplot.SkyLine` \
                                                     or :class:`~.gbplot.SkyPolygon`):
            The localization contour plot elements
        loc_posterior (:class:`~.gbmplot.SkyHeatmap`):
            The localization gradient plot element
        sun (:class:`~.gbmplot.Sun`): The Sun plot element
        text_color (str): The color of the text labels
    """
    _background = 'antiquewhite'
    _textcolor = 'black'
    _canvascolor = 'white'
    _x_origin = 180
    _y_origin = 0
    _fontsize = 10

    def __init__(self, canvas=None, projection='mollweide', flipped=True,
                 **kwargs):
        super().__init__(figsize=(10, 5), canvas=canvas, projection=projection,
                         **kwargs)

        # set up the plot background color and the sky grid
        self._figure.set_facecolor(self._canvascolor)
        self._ax.set_facecolor(self._background)
        self._ax.grid(True, linewidth=0.5)

        # create the axes tick labels
        self._longitude_axis(flipped)
        self._latitude_axis()

        self._flipped = flipped
        self._sun = None
        self._earth = None
        self._detectors = Collection()
        self._galactic_plane = None
        self._posterior = None
        self._clevels = Collection()

    @property
    def sun(self):
        return self._sun

    @property
    def earth(self):
        return self._earth

    @property
    def detectors(self):
        return self._detectors

    @property
    def gc(self):
        return self._galactic_plane

    @property
    def loc_posterior(self):
        return self._posterior

    @property
    def loc_contours(self):
        return self._clevels

    @property
    def canvas_color(self):
        return self._canvascolor

    @canvas_color.setter
    def canvas_color(self, color):
        self._figure.set_facecolor(color)
        self._canvascolor = color

    @property
    def background_color(self):
        return self._background

    @background_color.setter
    def background_color(self, color):
        self.ax.set_facecolor(color)
        self._background = color

    @property
    def text_color(self):
        return self._textcolor

    @text_color.setter
    def text_color(self, color):
        self._ax.set_yticklabels(self._ytick_labels, fontsize=self._fontsize,
                                 color=color)
        self._ax.set_xticklabels(self._xtick_labels, fontsize=self._fontsize,
                                 color=color)
        self._textcolor = color

    @property
    def fontsize(self):
        return self._fontsize

    @fontsize.setter
    def fontsize(self, size):
        self._ax.set_yticklabels(self._ytick_labels, fontsize=size,
                                 color=self._textcolor)
        self._ax.set_xticklabels(self._xtick_labels, fontsize=size,
                                 color=self._textcolor)
        self._fontsize = size

    def add_poshist(self, data, trigtime=None, detectors='all', geo=True,
                    sun=True,
                    galactic_plane=True):
        """Add a Position History or Trigdat object to plot the location of the
        Earth, Sun, and detector pointings

        Args:
            data (:class:`~gbm.data.PosHist` or :class:`~gbm.data.Trigdat`)
                A Position History or Trigdat object
            trigtime (float, optional): If data is PosHist, set trigtime to a 
                                        particular time of interest.
                                        The Trigdat trigger time overrides this
            detectors ('all' or list): A list of detectors or "all" to plot the 
                                   pointings on the sky
            geo (bool, optional): If True, plot the Earth. Default is True.
            sun (bool, optional): If True, plot the Sun. Default is True.
            galactic_plane (bool, optional):
                If True, plot the Galactic plane. Default is True.
        """
        if hasattr(data, 'trigtime'):
            trigtime = data.trigtime

        if detectors == 'all':
            dets = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6',
                    'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
        else:
            dets = detectors

        if trigtime is not None:
            if sun:
                sun_loc = get_sun_loc(trigtime)
                self.plot_sun(*sun_loc)
            if geo:
                geo_ra, geo_dec = data.get_geocenter_radec(trigtime)
                radius = data.get_earth_radius(trigtime)
                self.plot_earth(geo_ra, geo_dec, radius)
            for det in dets:
                ra, dec = data.detector_pointing(det, trigtime)
                self.plot_detector(ra, dec, det)

        # testing
        # lon = data.get_longitude(trigtime)
        # lat = data.get_latitude(trigtime)
        # alt = data.get_altitude(trigtime)
        # test_footprint(self.ax, self._earth._artists[0], lon, lat, alt,
        #               geo_ra, geo_dec, radius)

        if galactic_plane:
            self.plot_galactic_plane()

    def add_healpix(self, hpx, gradient=True, clevels=None, sun=True,
                    earth=True,
                    detectors='all', galactic_plane=True):
        """Add HealPix object to plot a localization and optionally the location 
        of the Earth, Sun, and detector pointings

        Args:
            hpx (:class:`~gbm.data.HealPix`): The HealPix object
            gradient (bool, optional): 
                If True, plots the posterior as a color gradient. If False, 
                plot the posterior as color-filled confidence regions.
            clevels (list of float, optional):
                The confidence levels to plot contours. By default plots at
                the 1, 2, and 3 sigma level.
            detectors ('all' or list):
                A list of detectors or "all" to plot the pointings on the sky
            earth (bool, optional): If True, plot the Earth. Default is True.
            sun (bool, optional): If True, plot the Sun. Default is True.
            galactic_plane (bool, optional):
                If True, plot the Galactic plane. Default is True.
        
        Note:
            Setting `gradient=False` when plotting an annulus may produce 
            unexpected results at this time.  It is suggested to use 
            `gradient=True` for plotting annuli maps.
        """
        if clevels is None:
            clevels = [0.997, 0.955, 0.687]
        if detectors == 'all':
            detectors = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6',
                         'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
        
        # determine what the resolution of the sky grid should be based on the
        # resolution of the healpix
        approx_res = np.sqrt(hpx.pixel_area)
        numpts_ra = int(np.floor(0.5*360.0/approx_res))
        numpts_dec = int(np.floor(0.5*180.0/approx_res))
        
        if gradient:
            prob_arr, ra_arr, dec_arr = hpx.prob_array(numpts_ra=numpts_ra,
                                                       numpts_dec=numpts_dec)
            self._posterior = self.plot_heatmap(prob_arr, ra_arr, dec_arr)

        for clevel in clevels:
            paths = hpx.confidence_region_path(clevel, numpts_ra=numpts_ra, 
                                               numpts_dec=numpts_dec)
            numpaths = len(paths)
            if gradient:
                for i in range(numpaths):
                    contour = SkyLine(paths[i][:, 0], paths[i][:, 1], self.ax,
                                      color='black',
                                      alpha=0.7, linewidth=2,
                                      flipped=self._flipped)
                    self._clevels.insert(str(clevel) + '_' + str(i), contour)
            else:
                for i in range(numpaths):
                    contour = SkyPolygon(paths[i][:, 0], paths[i][:, 1],
                                         self.ax, color='purple',
                                         face_alpha=0.3, flipped=self._flipped)
                    self._clevels.insert(str(clevel) + '_' + str(i), contour)

        # plot sun
        if sun:
            try:
                self.plot_sun(*hpx.sun_location)
            except:
                pass
        
        # plot earth
        if earth:
            try:
                geo_rad = 67.0 if hpx.geo_radius is None else hpx.geo_radius
                self.plot_earth(*hpx.geo_location, geo_rad)
            except:
                pass
        
        # plot detector pointings
        try:
            for det in detectors:
                self.plot_detector(*getattr(hpx, det + '_pointing'), det)
        except:
            pass
        
        # plot galactic plane
        if galactic_plane:
            self.plot_galactic_plane()

    def plot_sun(self, x, y, **kwargs):
        """Plot the sun

        Args:
            x (float): The RA of the Sun
            y (float): The Dec of the Sun
            **kwargs: Options to pass to :class:`~.gbmplot.Sun`
        """
        self._sun = Sun(x, y, self.ax, flipped=self._flipped, **kwargs)

    def plot_earth(self, x, y, radius, **kwargs):
        """Plot the Earth

        Args:
            x (float): The RA of the geocenter
            y (float): The Dec of the geocenter
            radius (float): The radius of the Earth, in degrees
            **kwargs: Options to pass to :class:`~.gbmplot.SkyCircle`
        """
        self._earth = SkyCircle(x, y, radius, self.ax, flipped=self._flipped,
                                color='deepskyblue', face_alpha=0.25,
                                edge_alpha=0.50, **kwargs)

    def plot_detector(self, x, y, det, radius=10.0, **kwargs):
        """Plot a detector pointing

        Args:
            x (float): The RA of the detector normal
            y (float): The Dec of the detector normal
            det (str): The detector name
            radius (float, optional): The radius of pointing, in degrees. 
                                      Default is 10.0
            **kwargs: Options to pass to :class:`~.gbmplot.SkyCircle`
        """
        pointing = DetectorPointing(x, y, radius, det, self.ax,
                                    flipped=self._flipped, **kwargs)
        self._detectors.insert(det, pointing)

    def plot_galactic_plane(self):
        """Plot the Galactic plane
        """
        self._galactic_plane = GalacticPlane(self.ax, flipped=self._flipped)

    def plot_heatmap(self, heatmap, ra_array, dec_array, **kwargs):
        """Plot a heatmap on the sky

        Args:
            heatmap (np.array): A 2D array of values
            ra_array (np.array): The array of RA gridpoints
            dec_array (np.array): The array of Dec gridpoints
            radius (float): The radius of pointing, in degrees
            **kwargs: Options to pass to :class:`~.gbmplot.SkyHeatmap`
        """
        heatmap = SkyHeatmap(ra_array, dec_array, heatmap, self.ax,
                             flipped=self._flipped, **kwargs)
        return heatmap

    def _longitude_axis(self, flipped):
        # longitude labels
        # these have to be shifted on the plot because matplotlib natively
        # goes from -180 to +180
        tick_labels = np.array(
            [210, 240, 270, 300, 330, 0, 30, 60, 90, 120, 150])
        tick_labels = tick_labels - 360 - int(self._x_origin)
        # flip coordinates
        if flipped:
            tick_labels = -tick_labels
        tick_labels = np.remainder(tick_labels, 360)

        # format the tick labels with degrees
        self._xtick_labels = [str(t) + '$^\circ$' for t in tick_labels]
        self._ax.set_xticklabels(self._xtick_labels, fontsize=self._fontsize,
                                 color=self._textcolor)

    def _latitude_axis(self):
        # latitude labels
        # matplotlib natively plots from -90 to +90 from bottom to top. 
        # this is fine for equatorial coordinates, but we have to shift if 
        # we are plotting in spacecraft coordinates
        tick_labels = np.array(
            [75, 60, 45, 30, 15, 0, -15, -30, -45, -60, -75])
        tick_labels = (self._y_origin - tick_labels)
        if np.sign(self._y_origin) == -1:
            tick_labels *= -1
        self._ytick_labels = [str(t) + '$^\circ$' for t in tick_labels]
        self._ax.set_yticklabels(self._ytick_labels, fontsize=self._fontsize,
                                 color=self._textcolor)


class FermiSkyPlot(SkyPlot):
    """Class for plotting in Fermi spacecraft coordinates.

    Parameters:
        projection (str, optional): The projection of the map. 
                                    Default is 'mollweide'
        **kwargs: Options to pass to :class:`~.gbmplot.GbmPlot`

    Attributes:
        ax (:class:`matplotlib.axes`): The matplotlib axes object for the plot
        background_color (str): The color of the plot background. This attribute
                                can be set.
        canvas (Canvas Backend object): The plotting canvas, if set upon 
                                        initialization.
        canvas_color (str): The color of the plotting canvas. This attribute 
                            can be set.
        detectors (:class:`~.gbmplot.Collection` of :class:`~.gbmplot.SkyCircle`):
            The collection of detector plot elements
        earth (:class:`~.gbmplot.SkyCircle`): The Earth plot element
        fig (:class:`matplotlib.figure`): The matplotlib figure object
        fontsize (int): The font size of the text labels. This attribute can be set.
        gc (:class:`~.gbmplot.GalacticPlane`):
            The reference to the galactic plane plot element
        loc_contours: (:class:`~.gbmplot.Collection` of :class:`~.gbmplot.SkyLine` \
                                                     or :class:`~.gbplot.SkyPolygon`):
            The localization contour plot elements
        loc_posterior (:class:`~.gbmplot.SkyHeatmap`):
            The localization gradient plot element
        sun (:class:`~.gbmplot.Sun`): The Sun plot element
        text_color (str): The color of the text labels
    """
    _y_origin = -90
    _x_origin = 0

    def __init__(self, canvas=None, projection='mollweide', **kwargs):
        super(FermiSkyPlot, self).__init__(canvas=canvas, flipped=False,
                                           projection=projection, **kwargs)

    def add_poshist(self, data, trigtime=None, detectors='all', geo=True,
                    sun=True, galactic_plane=True):
        if hasattr(data, 'trigtime'):
            trigtime = data.trigtime

        if detectors == 'all':
            dets = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6',
                    'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']

        if trigtime is not None:
            if sun:
                sun_loc = get_sun_loc(trigtime)
                sun_loc = data.to_fermi_frame(*sun_loc, trigtime)
                self.plot_sun(*sun_loc)
            if geo:
                ra, dec = data.get_geocenter_radec(trigtime)
                az, zen = data.to_fermi_frame(ra, dec, trigtime)
                radius = data.get_earth_radius(trigtime)
                self.plot_earth(az, zen, radius)

            for det in dets:
                self.plot_detector(det)

            if galactic_plane:
                quat = data.get_quaternions(trigtime)
                self.plot_galactic_plane(quat)

    def plot_detector(self, det, radius=10.0, **kwargs):
        """Plot a detector pointing

        Args:
            det (str): The detector name
            radius (float, optional): The radius of pointing, in degrees. 
                                      Default is 10.0
            **kwargs: Options to pass to :class:`~.gbmplot.SkyCircle`
        """
        super(FermiSkyPlot, self).plot_detector(0.0, 0.0, det, radius=radius,
                                                fermi=True, **kwargs)

    def plot_earth(self, x, y, radius, **kwargs):
        super(FermiSkyPlot, self).plot_earth(x, y, radius, fermi=True,
                                             **kwargs)

    def plot_sun(self, x, y, **kwargs):
        super(FermiSkyPlot, self).plot_sun(x, y, fermi=True, **kwargs)

        # mark TODO: enable loc plotting

    # will need the quaternion
    def add_healpix(*args, **kwargs):
        """Not yet implemented"""
        raise NotImplementedError('Not yet implemented for the Fermi frame')

    def plot_galactic_plane(self, quat):
        """Plot the Galactic plane
        
        Args:
            quat (np.array): The Fermi attitude quaternion
        """
        self._galactic_plane = GalacticPlane(self.ax, flipped=self._flipped,
                                             fermi_quat=quat)
