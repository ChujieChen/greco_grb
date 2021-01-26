# coords.py: Module containing coordinate conversion functions
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

import astropy.coordinates as coordinates
import astropy.time
import numpy as np
from astropy.utils.data import download_file
from astropy.utils.iers import IERS_A, IERS_A_URL
from .time import TimeFermiSec

def haversine(lon1, lat1, lon2, lat2, deg=True):
    """Calculates the angular separation between two points using the
    haversine equation. If degrees are passed, degrees are returned. else
    the input/output is assumed to be radians.
    lon -> azimuth
    lat -> zenith
    
    Args:
        lon1 (float): lon/az of first point
        lat1 (float): lat/zen of first point
        lon2 (float): lon/az of second point
        lat2 (float): lat/zen of second point
        deg (bool, optional): True if input/output in degrees.
    
    Returns:
        float: Angular separation between points
    """
    if deg:
        lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    d_lat = 0.5 * (lat2 - lat1)
    d_lon = 0.5 * (lon2 - lon1)

    a = np.sin(d_lat) ** 2 + (np.sin(d_lon) ** 2 * np.cos(lat1) * np.cos(lat2))
    alpha = 2. * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    if deg:
        alpha = np.rad2deg(alpha)

    return alpha


def azzen_to_cartesian(az, zen, deg=True):
    """Convert spacecraft azimuth/zenith to Cartesian coordinates
    
    Args:
        az (float or np.array): Spacecraft azimuth
        zen (float or np.array): Spacecraft zenith
        deg (bool, optional): True (default) if input is in degrees, 
                              otherwise input is radians
    
    Returns:
        np.array: 3-element Cartesian coordinate
    """
    if deg:
        az = np.deg2rad(az)
        zen = np.deg2rad(zen)

    sinZen = np.sin(zen)
        
    return np.array(
        [sinZen * np.cos(az), sinZen * np.sin(az), np.cos(zen)])


def radec_to_cartesian(ra, dec, deg=True):
    """Convert RA/Dec to Cartesian coordinates
    
    Args:
        ra (float or np.array): Right Ascension
        dec (float or np.array): Declination
        deg (bool, optional): True (default) if input is in degrees, 
                              otherwise input is radians
    
    Returns:
        np.array: 3-element Cartesian coordinate
    """
    if deg:
        ra = np.deg2rad(ra)
        dec = np.deg2rad(dec)

    return np.array(
        [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def quaternion_conj(quat):
    """Calculate conjugate of quaternion
    
    Note: 
        GBM quaternions are defined with the last element as the scalar value
    
    Args:
        quat (np.array): 4-element quaternion
    
    Returns:
        np.array: quaternion conjugate
    """
    cquat = np.copy(quat)
    cquat[0:3] = -cquat[0:3]
    return cquat


def quaternion_inv(quat):
    """Calculate inverse of quaternion
    
    Note: 
        GBM quaternions are defined with the last element as the scalar value

    Args:
        quat (np.array): 4-element quaternion
    
    Returns:
        np.array: quaternion inverse
    """
    cquat = quaternion_conj(quat)
    iquat = cquat / np.sum(quat ** 2)
    return iquat


def quaternion_prod(quat1, quat2):
    """Calculate product of two quaternions
    
    Note:
        GBM quaternions are defined with the last element as the scalar value

    Args:
        quat1 (np.array): 4-element quaternion
        quat2 (np.array): 4-element quaternion
    
    Returns:
        np.array: product of the quaternions
    """
    q = np.copy(quat1)
    r = np.copy(quat2)
    q[0] = quat1[3]
    q[1:] = quat1[0:3]
    r[0] = quat2[3]
    r[1:] = quat2[0:3]

    t = np.copy(quat1)
    t[0] = r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3]
    t[1] = r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2]
    t[2] = r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1]
    t[3] = r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0]
    quatprod = np.copy(t)
    quatprod[3] = t[0]
    quatprod[0:3] = t[1:]
    return quatprod


def spacecraft_direction_cosines(quat):
    """Convert `n` spacecraft quaternions to direction cosine matrix
    
    Args:
        quat (np.array): (4, `n`) element array of quaternions
        
    Returns:
        np.array: (3,3, `n`) array containing the spacecraft direction cosines
    """
    n_elements = quat.shape[0]
    if n_elements != 4:
        raise ValueError(
            'Quaternion must have 4 elements, {0} given'.format(n_elements))

    ndim = len(quat.shape)
    if ndim == 2:
        numquats = quat.shape[1]
    else:
        numquats = 1
    quat /= np.linalg.norm(quat, axis=0)

    sc_cosines = np.zeros((3, 3, numquats), dtype=float)
    sc_cosines[0, 0, np.newaxis] = (
                quat[0, np.newaxis] ** 2 - quat[1, np.newaxis] ** 2 - quat[
            2, np.newaxis] ** 2 +
                quat[3, np.newaxis] ** 2)
    sc_cosines[1, 0, np.newaxis] = 2.0 * (
                quat[0, np.newaxis] * quat[1, np.newaxis] + quat[
            3, np.newaxis] *
                quat[2, np.newaxis])
    sc_cosines[2, 0, np.newaxis] = 2.0 * (
                quat[0, np.newaxis] * quat[2, np.newaxis] - quat[
            3, np.newaxis] *
                quat[1, np.newaxis])
    sc_cosines[0, 1, np.newaxis] = 2.0 * (
                quat[0, np.newaxis] * quat[1, np.newaxis] - quat[
            3, np.newaxis] *
                quat[2, np.newaxis])
    sc_cosines[1, 1, np.newaxis] = (
                -quat[0, np.newaxis] ** 2 + quat[1, np.newaxis] ** 2 - quat[
            2, np.newaxis] ** 2 +
                quat[3, np.newaxis] ** 2)
    sc_cosines[2, 1, np.newaxis] = 2.0 * (
                quat[1, np.newaxis] * quat[2, np.newaxis] + quat[
            3, np.newaxis] *
                quat[0, np.newaxis])
    sc_cosines[0, 2, np.newaxis] = 2.0 * (
                quat[0, np.newaxis] * quat[2, np.newaxis] + quat[
            3, np.newaxis] *
                quat[1, np.newaxis])
    sc_cosines[1, 2, np.newaxis] = 2.0 * (
                quat[1, np.newaxis] * quat[2, np.newaxis] - quat[
            3, np.newaxis] *
                quat[0, np.newaxis])
    sc_cosines[2, 2, np.newaxis] = (
                -quat[0, np.newaxis] ** 2 - quat[1, np.newaxis] ** 2 + quat[
            2, np.newaxis] ** 2 +
                quat[3, np.newaxis] ** 2)
    return np.squeeze(sc_cosines)


def geocenter_direction_cosines(gz, az, zen, deg=True):
    """Create `n` geocenter direction cosine matrix(es) for a position or 
    positions in spacecraft coordinates
    
    Args:
        gz (np.array): Geocenter vector in spacecraft Cartesian coordinates
        az (float or np.array): The azimuth location of positions in spacecraft 
                                coordinates
        zen (float or np.array): The zenith location of positions in spacecraft 
                                 coordinates    
        deg (bool, optional): If True, the input is in degrees. Default is True.
    
    Returns:
        np.array: (3,3, `n`) array containing the geocenter direction cosines \
                  for each point
    """
    try:
        numpts = len(az)
    except:
        az = np.asarray([az])
        zen = np.asarray([zen])
        numpts = len(az)

    geo_cosines = np.zeros((3, 3, numpts), dtype=float)

    # the vector from the source position to the geocenter defines our frame.
    # specifically, the source->geocenter vector is our z-axis

    # geocenter vector in cartesian sc coordinates
    # gz = azzen_to_cartesian(geo_az, geo_zen)
    gz /= np.linalg.norm(gz)  # np.sqrt(np.sum(gz**2))

    # source vector in cartesian sc coordinates
    sl = azzen_to_cartesian(az, zen)
    sl /= np.linalg.norm(sl, axis=0)

    # y-axis of our frame new frame
    gy = np.array([gz[1] * sl[2, :] - gz[2] * sl[1, :],
                   gz[2] * sl[0, :] - gz[0] * sl[2, :],
                   gz[0] * sl[1, :] - gz[1] * sl[0, :]])
    gy /= np.linalg.norm(gy, axis=0)

    # x-axis of our new frame
    gx = np.array([gy[1, :] * gz[2] - gy[2, :] * gz[1],
                   gy[2, :] * gz[0] - gy[0, :] * gz[2],
                   gy[0, :] * gz[1] - gy[1, :] * gz[0]])
    gx /= np.linalg.norm(gx, axis=0)

    # complete direction cosine matrix
    geo_cosines[0, :, :] = gx
    geo_cosines[1, :, :] = gy
    geo_cosines[2, :, :] = gz[:, np.newaxis]
    return np.squeeze(geo_cosines)


def geocenter_to_spacecraft(phi, theta, geo_cosines):
    """Convert `n` points from geocenter frame to spacecraft frame
    
    Args:
        phi (float or np.array): The azimuthal location of positions in 
                                 geocentric coordinates
        theta (float or np.array): The polar location of positions in 
                                   geocentric coordinates    
        geo_cosines (np.array): (3,3) array containing the direction cosines
    
    Returns:
        (np.array, np.array): Azimuth and zenith positions in the spacecraft frame
    """
    try:
        numpts = len(phi)
    except:
        phi = np.asarray([phi])
        theta = np.asarray([theta])
        numpts = len(phi)

    # spherical coordinate vector in cartesian coordinates
    a = azzen_to_cartesian(phi, theta)

    gx = geo_cosines[0, :]
    gy = geo_cosines[1, :]
    gz = geo_cosines[2, :]

    # the cartesian vector rotated into the spacecraft frame
    dir = [a[0, :] * gx[0] + a[1, :] * gy[0] + a[2, :] * gz[0],
           a[0, :] * gx[1] + a[1, :] * gy[1] + a[2, :] * gz[1],
           a[0, :] * gx[2] + a[1, :] * gy[2] + a[2, :] * gz[2]]
    dir = -np.array(dir)
    dir /= np.linalg.norm(dir, axis=0)

    # convert cartesian vector to az/zen
    az = np.rad2deg(np.arctan2(dir[1, :], dir[0, :]))
    az[az < 0.0] += 360.0
    zen = np.rad2deg(np.arccos(dir[2, :]))
    return (np.squeeze(az), np.squeeze(zen), np.squeeze(dir))


def geocenter_in_radec(coord):
    """Calculate the location of the Earth center RA and Dec
    
    Args:
        coord (np.array): (3, `n`) array containing Geocentric cartesian 
                          coordinates in meters
    
    Returns:
        (np.array, np.array): RA and Dec of Earth center as viewed by \
                              Fermi in degrees
    """
    unit_vec = -coord / np.linalg.norm(-coord, axis=0)
    dec = np.pi / 2.0 - np.arccos(unit_vec[2, np.newaxis])
    ra = np.arctan2(unit_vec[1, np.newaxis], unit_vec[0, np.newaxis])
    ra[ra < 0.0] += 2.0 * np.pi
    return np.squeeze(np.rad2deg(ra)), np.squeeze(np.rad2deg(dec))


def spacecraft_to_radec(az, zen, quat, deg=True):
    """Convert a position in spacecraft coordinates (Az/Zen) to J2000 RA/Dec
    The options for input for this function are as follows:
    
    * a single Az/Zen position and multiple attitude transforms
    * multiple Az/Zen positions and a single attitude transform
    * multiple Az/Zen positions each with a corresponding attitude transform
    
    Args:
        az (float or np.array): Spacecraft azimuth  
        zen (float or np.array): Spacecraft zenith
        quat (np.array): (4, `n`) spacecraft attitude quaternion array
        deg (bool, optional): True if input/output in degrees.
    
    Returns:
        (np.array, np.array): RA and Dec of the transformed position
    """
    ndim = len(quat.shape)
    if ndim == 2:
        numquats = quat.shape[1]
    else:
        numquats = 1

    # convert az/zen to cartesian coordinates
    pos = azzen_to_cartesian(az, zen, deg=deg)
    ndim = len(pos.shape)
    if ndim == 2:
        numpos = pos.shape[1]
    else:
        numpos = 1

    # spacecraft direction cosine matrix
    sc_cosines = spacecraft_direction_cosines(quat)

    # can do one sky position over many transforms, many sky positions over one
    # transform, or a transform for each sky position
    if (numpos == 1) & (numquats > 1):
        pos = np.repeat(pos, numquats).reshape(3, -1)
        numdo = numquats
    elif (numpos > 1) & (numquats == 1):
        sc_cosines = np.repeat(sc_cosines, numpos).reshape(3, 3, -1)
        numdo = numpos
    elif numpos == numquats:
        numdo = numpos
        if numdo == 1:
            sc_cosines = sc_cosines[:, :, np.newaxis]
            pos = pos[:, np.newaxis]
    else:
        raise ValueError(
            'If the size of az/zen coordinates is > 1 AND the sizeof quaternions is > 1, then they must be of same size'
        )

    # convert numpy arrays to list of arrays for vectorized calculations
    sc_cosines_list = np.squeeze(np.split(sc_cosines, numdo, axis=2))
    pos_list = np.squeeze(np.split(pos, numdo, axis=1))
    if numdo == 1:
        sc_cosines_list = [sc_cosines_list]
        pos_list = [pos_list]

    # convert position to J2000 frame
    cartesian_pos = np.array(list(map(np.dot, sc_cosines_list, pos_list))).T
    cartesian_pos[2, (cartesian_pos[2, np.newaxis] < -1.0).reshape(-1)] = -1.0
    cartesian_pos[2, (cartesian_pos[2, np.newaxis] > 1.0).reshape(-1)] = 1.0

    # transform cartesian position to RA/Dec in J2000 frame
    dec = np.arcsin(cartesian_pos[2, np.newaxis])
    ra = np.arctan2(cartesian_pos[1, np.newaxis], cartesian_pos[0, np.newaxis])
    ra[(np.abs(cartesian_pos[1, np.newaxis]) < 1e-6) & (
                np.abs(cartesian_pos[0, np.newaxis]) < 1e-6)] = 0.0
    ra[ra < 0.0] += 2.0 * np.pi

    if deg:
        ra = np.rad2deg(ra)
        dec = np.rad2deg(dec)
    
    return np.squeeze(ra), np.squeeze(dec)


def radec_to_spacecraft(ra, dec, quat, deg=True):
    """Convert a position in J2000 RA/Dec to spacecraft coordinates (Az/Zen).
    The options for input for this function are as follows:
    
    * a single RA/Dec position and multiple attitude transforms
    * multiple RA/Dec positions and a single attitude transform
    * multiple RA/Dec positions each with a corresponding attitude transform
    
    Args:
        ra (float or np.array): Right Ascension  
        dec (float or np.array): Declination
        quat (np.array): (4, `n`) spacecraft attitude quaternion array
        deg (bool, optional): True if input/output in degrees.

    Returns:
        (np.array, np.array): Spacecraft azimuth and zenith of the transformed \
                              position
    """
    ndim = len(quat.shape)
    if ndim == 2:
        numquats = quat.shape[1]
    else:
        numquats = 1

    # convert az/zen to cartesian coordinates
    pos = radec_to_cartesian(ra, dec, deg=deg)
    ndim = len(pos.shape)
    if ndim == 2:
        numpos = pos.shape[1]
    else:
        numpos = 1

    # spacecraft direction cosine matrix
    sc_cosines = spacecraft_direction_cosines(quat)

    # can do one sky position over many transforms, many sky positions over one
    # transform, or a transform for each sky position
    if (numpos == 1) & (numquats > 1):
        pos = np.repeat(pos, numquats).reshape(3, -1)
        numdo = numquats
    elif (numpos > 1) & (numquats == 1):
        sc_cosines = np.repeat(sc_cosines, numpos).reshape(3, 3, -1)
        numdo = numpos
    elif numpos == numquats:
        numdo = numpos
        if numdo == 1:
            sc_cosines = sc_cosines[:, :, np.newaxis]
            pos = pos[:, np.newaxis]
    else:
        raise ValueError(
            'If the size of az/zen coordinates is > 1 AND the sizeof quaternions is > 1, then they must be of same size'
        )

    # convert numpy arrays to list of arrays for vectorized calculations
    sc_cosines_list = np.squeeze(np.split(sc_cosines, numdo, axis=2))
    pos_list = np.squeeze(np.split(pos, numdo, axis=1))
    if numdo == 1:
        sc_cosines_list = [sc_cosines_list]
        pos_list = [pos_list]

    # convert position from J2000 frame    
    cartesian_pos = np.array(list(map(np.dot, pos_list, sc_cosines_list))).T

    # convert Cartesian coordinates to spherical
    zen = np.arccos(cartesian_pos[2, np.newaxis])
    az = np.arctan2(cartesian_pos[1, np.newaxis], cartesian_pos[0, np.newaxis])
    az[(np.abs(cartesian_pos[1, np.newaxis]) < 1e-6) & (
                np.abs(cartesian_pos[0, np.newaxis]) < 1e-6)] = 0.0
    az[az < 0.0] += 2.0 * np.pi

    if deg:
        az = np.rad2deg(az)
        zen = np.rad2deg(zen)
    
    return np.squeeze(az), np.squeeze(zen)


def latitude_from_geocentric_coords_simple(coord):
    """Calculate latitude from Geocentric Cartesian coordinates. Assumes the 
    Earth is a simple sphere.  Also returns altitude.
    
    Args:
        coord (np.array): (3, `n`) array containing Geocentric cartesian 
                          coordinates in meters
    
    Returns:
        (np.array, np.array): Position of Fermi in Earth Latitude and altitude \
                              in meters
    """
    mean_earth_radius = 6371009.0  # m
    radius = np.sqrt(np.sum(coord[:, np.newaxis] ** 2, axis=0))
    latitude = np.rad2deg(np.arcsin(coord[2, np.newaxis] / radius))
    altitude = radius - mean_earth_radius
    return np.squeeze(latitude), np.squeeze(altitude)


def latitude_from_geocentric_coords_complex(coord):
    """Calculate latitude from Geocentric Cartesian coordinates. Uses the 
    WGS 1984 model of the shape of the Earth.  Also returns altitude.
    
    Args:
        coord (np.array): (3, `n`) array containing Geocentric cartesian 
                          coordinates in meters
    
    Returns:
        (np.array, np.array): Position of Fermi in Earth Latitude and altitude \
                              in meters
    """
    ndim = len(coord.shape)
    if ndim == 2:
        numpoints = coord.shape[1]
    else:
        numpoints = 1

    # Parameters of the World Geodetic System 1984
    # semi-major axis
    wgs84_a = 6378137.0  # km
    # reciprocal of flattening
    wgs84_1overf = 298.257223563

    rho = np.sqrt(coord[0, np.newaxis] ** 2 + coord[1, np.newaxis] ** 2)
    f = 1.0 / wgs84_1overf
    e_sq = 2.0 * f - f ** 2

    # should completely converge in 3 iterations
    n_iter = 3
    kappa = np.zeros((n_iter + 1, numpoints), dtype=float)
    kappa[0, :np.newaxis] = 1.0 / (1.0 - e_sq)
    for i in range(1, n_iter + 1):
        c = (rho ** 2 + (1.0 - e_sq) * coord[2, np.newaxis] ** 2 * kappa[
            i - 1, np.newaxis] ** 2) ** 1.5 / (
                    wgs84_a * e_sq)
        kappa[i, np.newaxis] = (c + (1.0 - e_sq) * coord[2, np.newaxis] ** 2 *
                                kappa[i - 1, np.newaxis] ** 3) / (
                                       c - rho ** 2)

    phi = np.arctan(kappa[-1, np.newaxis] * coord[2, np.newaxis] / rho)
    h = (1.0 / kappa[-1, np.newaxis] - 1.0 / kappa[0, np.newaxis]) * \
        np.sqrt(rho ** 2 + coord[2, np.newaxis] ** 2 * kappa[
            -1, np.newaxis] ** 2) / e_sq
    latitude = np.rad2deg(phi)
    altitude = h
    return np.squeeze(latitude), np.squeeze(altitude)


def longitude_from_geocentric_coords(coord, met, ut1=False):
    """Calculate the East longitude from Geocentric coordinates. Requires the 
    conversion of Fermi MET to sidereal time to rotate the Earth under the 
    spacecraft.  The conversion can either be performed using UTC (less 
    accurate) or UT1 (more accurate) which uses IERS tables to correct for 
    variations in Earth rotation velocity
    
    Args:
        coord (np.array): (3, `n`) array containing Geocentric cartesian 
                          coordinates in meters
        met: (float or np.array): The MET time(s)
        ut1 (bool, optional): If True, use UT1 instead of UTC to calculate 
                              sidereal time. Default is False
    
    Returns:
        np.array or float: Position of Fermi in East Longitude
    """
    ndim = len(coord.shape)
    if ndim == 2:
        if coord.shape[1] != met.shape[0]:
            raise ValueError(
                'The number of coordinate positions must match the number of times')

    # get Fermi time object in UTC standard; need astropy for sidereal time conversion
    utc_time = astropy.time.Time(met, format='fermi').utc
    if ut1 is True:
        try:
            # Use IERS table for sidereal time conversion
            utc_time.delta_ut1_utc = utc_time.get_delta_ut1_utc()
            theta_deg = utc_time.sidereal_time('apparent', 'greenwich').degree
        except:
            # IERS table is not current.  Update IERS table
            print("IERS table is not current.  Checking for updated Table.")
            iers_a_file = download_file(IERS_A_URL, cache=True)
            iers_a = IERS_A.open(iers_a_file)
            utc_time.delta_ut1_utc = utc_time.get_delta_ut1_utc(iers_a)
            theta_deg = utc_time.sidereal_time('apparent', 'greenwich').degree
    else:
        # Use less accurate UTC for sidereal time conversion
        utc_time.delta_ut1_utc = 0.0
        theta_deg = utc_time.sidereal_time('apparent', 'greenwich').degree

    inertial_longitude = np.arctan2(coord[1, np.newaxis], coord[0, np.newaxis])
    east_longitude = np.rad2deg(inertial_longitude) - theta_deg

    east_longitude[east_longitude < 0.0] += 360.0
    east_longitude[east_longitude < 0.0] += 360.0
    return np.squeeze(east_longitude)


def calc_mcilwain_l(latitude, longitude):
    """Estimate the McIlwain L value given the latitude (-30, +30) and 
    East Longitude.  This uses a cubic polynomial approximation to the full 
    calculation and is similar to the approach used by the GBM FSW.
   
    Args:
        latitude (np.array): Latitude in degrees from -180 to 180
        longitude (np.array): East longitude in degrees from 0 to 360
    
    Returns:
        np.array: McIlwain L value
    """

    latitude = np.asarray([latitude])
    longitude = np.asarray([longitude])
    orig_shape = latitude.shape
    latitude = latitude.flatten()
    longitude = longitude.flatten()
    # numPts = latitude.shape[0]
    coeffs_file = os.path.join(os.path.dirname(__file__),
                               'McIlwainL_Coeffs.npy')
    poly_coeffs = np.load(coeffs_file)
    longitude[longitude < 0.0] += 360.0
    longitude[longitude == 360.0] = 0.0

    bad_idx = (latitude < -30.0) | (latitude > 30.0) | (longitude < 0.0) | (
                longitude >= 360.0)
    if np.sum(bad_idx) != 0:
        raise ValueError(
            'Out of range coordinates for McIlwain L for {0} locations'.format(
                np.sum(bad_idx)))

    idx = np.asarray((longitude / 10.0).astype(int))
    idx2 = np.asarray(idx + 1)
    idx2[idx2 >= 36] = 0
    idx2 = idx2.astype(int)

    longitude_left = 10.0 * idx
    f = (longitude - longitude_left) / 10.0  # interpolation weight, 0 to 1

    try:
        num_pts = len(latitude)
    except:
        num_pts = 1
    mc_l = np.zeros(num_pts)
    for i in range(num_pts):
        mc_l[i] = (1.0 - f[i]) * (
                poly_coeffs[idx[i], 0] + poly_coeffs[idx[i], 1] * latitude[i] +
                poly_coeffs[idx[i], 2] *
                latitude[i] ** 2 + poly_coeffs[idx[i], 3] * latitude[i] ** 3) + \
                  f[i] * (
                          poly_coeffs[idx2[i], 0] + poly_coeffs[idx2[i], 1] *
                          latitude[i] + poly_coeffs[idx2[i], 2] *
                          latitude[i] ** 2 + poly_coeffs[idx2[i], 3] *
                          latitude[i] ** 3)
    mc_l = mc_l.reshape(orig_shape)
    return np.squeeze(mc_l)


def get_sun_loc(met):
    """Calculate sun location in RA/Dec for a given MET.
    
    Args:
        met (float): The MET time(s)

    Returns:
        (float, float): RA and Dec of the sun
    """

    utc_time = astropy.time.Time(met, format='fermi').utc
    sun = coordinates.get_sun(utc_time)
    return sun.ra.degree, sun.dec.degree


def saa_boundary():
    """The coordinates of the SAA boundary in latitude and East longitude
    
    Returns:
        (np.array, np.array): The latitude and East longitude values
    """
    lat_saa = np.array([-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                        -6.155, -8.880, -14.220, -18.404, -30.000, -30.000])
    lon_saa = np.array(
        [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
         -84.000, -89.200, -94.300, -94.300, -86.100, 33.900])
    return (lat_saa, lon_saa)
