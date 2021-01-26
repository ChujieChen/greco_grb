# headers.py: GBM data header definitions
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
import astropy.io.fits as fits
from gbm.time import Met
from gbm.detectors import Detector
import gbm


class GBMDefinitions:
    def __init__(self):
        self.telescope = ('telescop', 'GLAST', 'Name of mission/satellite')
        self.instrument = (
        'instrume', 'GBM', 'Specific instrument used for observation')
        self.observer = ('observer', 'Meegan', 'GLAST Burst Monitor P.I.')
        self.origin = ('origin', 'GIOC', 'Name of organization making file')
        self.timesys = ('timesys', 'TT', 'Time system used in time keywords')
        self.timeunit = (
        'timeunit', 's', 'Time since MJDREF, used in TSTART and TSTOP')
        self.mjdrefi = (
        'mjdrefi', 51910, 'MJD of GLAST reference epoch, integer part')
        self.mjdreff = ('mjdreff', '7.428703703703703e-4',
                        'MJD of GLAST reference epoch, fractional part')
        self.radecsys = ('radecsys', 'FK5', 'Stellar reference frame')
        self.equinox = ('equinox', 2000.0, 'Equinox for RA and Dec')


def _set_str(val):
    return str(val) if val is not None else None


def _set_float(val):
    return float(val) if val is not None else None


def _set_int(val):
    return int(val) if val is not None else None


def _met_to_utc(met):
    if met is not None:
        _met = Met(met)
        return _met.iso()


def _creator(version):
    return ('creator',
            'GBM Data Tools {} Software and version creating file'.format(
                version))


def _filetype(ftype):
    return ('filetype', ftype, 'Name for this type of FITS file')


def _detnam(detnam):
    try:
        det = Detector.from_str(detnam).long_name
    except:
        det = detnam
    return ('detnam', det, 'Individual detector name')


def _date(date):
    return ('date', date, 'file creation date (YYYY-MM-DDThh:mm:ss UT)')


def _dateobs(date_obs):
    return ('date-obs', date_obs, 'Date of start of observation')


def _dateend(date_end):
    return ('date-end', date_end, 'Date of end of observation')


def _tstart(tstart):
    return ('tstart', tstart, '[GLAST MET] Observation start time')


def _tstop(tstop):
    return ('tstop', tstop, '[GLAST MET] Observation stop time')


def _trigtime(trigtime):
    return (
    'trigtime', trigtime, 'Trigger time relative to MJDREF, double precision')


def _object(object):
    return ('object', object, 'Burst name in standard format, yymmddfff')


def _raobj(ra_obj):
    return ('ra_obj', ra_obj, 'Calculated RA of burst')


def _decobj(dec_obj):
    return ('dec_obj', dec_obj, 'Calculated Dec of burst')


def _errrad(err_rad):
    return ('err_rad', err_rad, 'Calculated Location Error Radius')


def _infile01(infile):
    return ('infile01', infile, 'Level 1 input lookup table file')


def _extname(extname):
    return ('extname', extname, 'name of this binary table extension')


def _hduclass():
    return (
    'hduclass', 'OGIP', 'Conforms to OGIP standard indicated in HDUCLAS1')


def _hduvers():
    return ('hduvers', '1.2.1', 'Version of HDUCLAS1 format in use')


def _chantype():
    return ('chantype', 'PHA', 'No corrections have been applied')


def _filter():
    return ('filter', 'None', 'The instrument filter in use (if any)')


def _detchans(detchans):
    return ('detchans', detchans, 'Total number of channels in each rate')


def _extver():
    return ('extver', 1, 'Version of this extension format')


def _tzero_(i, trigtime):
    return ('tzero{}'.format(i), trigtime, 'Offset, equal to TRIGTIME')


def primary(detnam=None, filetype=None, tstart=None, tstop=None, filename=None,
            trigtime=None, object=None, ra_obj=None, dec_obj=None,
            err_rad=None,
            infile=None):
    # enforce strings
    detnam = _set_str(detnam)
    filetype = _set_str(filetype)
    filename = _set_str(filename)
    object = _set_str(object)
    infile = _set_str(infile)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_creator(gbm.__version__))
    header.append(_filetype(filetype))
    header.append(
        ('file-ver', '1.0.0', 'Version of the format for this filetype'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(('filename', filename, 'Name of this file'))
    header.append(_trigtime(trigtime))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))
    header.append(_infile01(infile))

    return header


def ebounds(detnam=None, tstart=None, tstop=None, trigtime=None, object=None,
            ra_obj=None, dec_obj=None, err_rad=None, detchans=None,
            ch2e_ver=None,
            gain=None, infile=None):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)
    ch2e_ver = _set_str(ch2e_ver)
    infile = _set_str(infile)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)
    gain = _set_float(gain)

    # enforce ints
    detchans = _set_int(detchans)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_extname('EBOUNDS'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_trigtime(trigtime))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))
    header.append(_hduclass())
    header.append(
        ('hduclas1', 'RESPONSE', 'These are typically found in RMF files'))
    header.append(('hduclas2', 'EBOUNDS', 'From CAL/GEN/92-002'))
    header.append(_hduvers())
    header.append(_chantype())
    header.append(_filter())
    header.append(_detchans(detchans))
    header.append(_extver())
    header.append(
        ('ch2e_ver', ch2e_ver, 'Channel to energy conversion scheme used'))
    header.append(
        ('gain_cor', gain, 'Gain correction factor applied to energy edges'))
    header.append(_infile01(infile))

    return header


def spectrum(detnam=None, tstart=None, tstop=None, trigtime=None, object=None,
             ra_obj=None, dec_obj=None, err_rad=None, detchans=None,
             poisserr=True):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)

    # enforce ints
    detchans = _set_int(detchans)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_tzero_(4, trigtime))
    header.append(_tzero_(5, trigtime))
    header.append(_extname('SPECTRUM'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_trigtime(trigtime))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))
    header.append(_filter())
    header.append(
        ('areascal', 1., 'No special scaling of effective area by channel'))
    header.append(
        ('backfile', 'none', 'Name of corresponding background file (if any)'))
    header.append(('backscal', 1., 'No scaling of background'))
    header.append(
        ('corrfile', 'none', 'Name of corresponding correction file (if any)'))
    header.append(('corrscal', 1., 'Correction scaling file'))
    header.append(
        ('respfile', 'none', 'Name of corresponding RMF file (if any)'))
    header.append(
        ('ancrfile', 'none', 'Name of corresponding ARF file (if any)'))
    header.append(('sys_err', 0., 'No systematic errors'))
    header.append(('poisserr', poisserr, 'Assume Poisson Errors'))
    header.append(('grouping', 0, 'No special grouping has been applied'))
    header.append(_hduclass())
    header.append(
        ('hduclas1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'))
    header.append(
        ('hduclas2', 'TOTAL', 'Indicates gross data (source + background)'))
    header.append(('hduclas3', 'COUNT', 'Indicates data stored as counts'))
    header.append(('hduclas4', 'TYPEII', 'Indicates PHA Type II file format'))
    header.append(_hduvers())
    header.append(_chantype())
    header.append(_detchans(detchans))
    header.append(_extver())

    return header


def events(detnam=None, tstart=None, tstop=None, trigtime=None, object=None,
           ra_obj=None, dec_obj=None, err_rad=None, detchans=None):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)

    # enforce ints
    detchans = _set_int(detchans)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_tzero_(1, trigtime))
    header.append(_extname('EVENTS'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_trigtime(trigtime))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))
    header.append(
        ('respfile', 'none', 'Name of corresponding RMF file (if any)'))
    header.append(('evt_dead', 2.6e-6, 'Deadtime per event (s)'))
    header.append(_detchans(detchans))
    header.append(_hduclass())
    header.append(('hduclas1', 'EVENTS', 'Extension contains Events'))
    header.append(_extver())

    return header


def gti(detnam=None, tstart=None, tstop=None, trigtime=None, object=None,
        ra_obj=None, dec_obj=None, err_rad=None):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_tzero_(1, trigtime))
    header.append(_tzero_(2, trigtime))
    header.append(_extname('GTI'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_hduclass())
    header.append(('hduclas1', 'GTI', 'Indicates good time intervals'))
    header.append(_hduvers())
    header.append(_extver())
    header.append(_trigtime(trigtime))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))

    return header


def specresp(detnam=None, tstart=None, tstop=None, trigtime=None, object=None,
             ra_obj=None, dec_obj=None, mat_type=None, rsp_num=None,
             src_az=None,
             src_el=None, geo_az=None, geo_el=None, det_ang=None, geo_ang=None,
             numebins=None, detchans=None, infiles=None, atscat=None):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)
    mat_type = _set_str(mat_type)
    if infiles is None:
        infiles = [None] * 3
    else:
        infiles = [str(infile) for infile in infiles]
    atscat = _set_str(atscat)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    src_az = _set_float(src_az)
    src_el = _set_float(src_el)
    geo_az = _set_float(geo_az)
    geo_el = _set_float(geo_el)
    det_ang = _set_float(det_ang)
    geo_ang = _set_float(geo_ang)

    # enforce ints
    detchans = _set_int(detchans)
    rsp_num = _set_int(rsp_num)
    numebins = _set_int(numebins)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_extname('SPECRESP MATRIX'))
    header.append(_extver())
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_trigtime(trigtime))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('mat_type', mat_type, 'Response Matrix Type'))
    header.append(('rsp_num', rsp_num, 'Response matrix index number'))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(
        ('src_az', src_az, 'Azimuth of source in spacecraft coordinates'))
    header.append(
        ('src_el', src_el, 'Elevation of source in spacecraft coordinates'))
    header.append(
        ('geo_az', geo_az, 'Azimuth of geocenter in spacecraft coordinates'))
    header.append(
        ('geo_el', geo_el, 'Elevation of geocenter in spacecraft coordinates'))
    header.append(
        ('det_ang', det_ang, 'Angle between source and detector normal'))
    header.append(
        ('geo_ang', geo_ang, 'Angle between geocenter and detector normal'))
    header.append(_filter())
    header.append(_chantype())
    header.append(
        ('numebins', numebins, 'Number of true energy bins of the MATRIX'))
    header.append(_detchans(detchans))
    header.append(('infile01', infiles[0], 'Detector response database in'))
    header.append(('infile02', infiles[1], 'Detector response database in'))
    header.append(('infile03', infiles[2], 'Detector response database in'))
    header.append(('infile04', atscat, 'Atmospheric scattering datab'))
    header.append(_hduclass())
    header.append(_hduvers())
    header.append(('hduclas1', 'RESPONSE', 'Typically found in RMF files'))
    header.append(('hduclas2', 'RSP_MATRIX', 'From CAL/GEN/92-002'))

    return header


def pha_spectrum(detnam=None, tstart=None, tstop=None, trigtime=None,
                 object=None,
                 ra_obj=None, dec_obj=None, err_rad=None, detchans=None,
                 poisserr=True,
                 datatype=None, backfile=None, rspfile=None, exposure=None):
    # enforce strings
    detnam = _set_str(detnam)
    object = _set_str(object)
    datatype = _set_str(datatype)
    backfile = _set_str(backfile)
    rspfile = _set_str(rspfile)

    # enforce floats
    tstart = _set_float(tstart)
    tstop = _set_float(tstop)
    trigtime = _set_float(trigtime)
    ra_obj = _set_float(ra_obj)
    dec_obj = _set_float(dec_obj)
    err_rad = _set_float(err_rad)
    exposure = _set_float(exposure)

    # enforce ints
    detchans = _set_int(detchans)

    # do MET -> UTC time conversion
    date_obs = _met_to_utc(tstart)
    date_end = _met_to_utc(tstop)

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    # header.append(_tzero_(4, trigtime))
    # header.append(_tzero_(5, trigtime))
    header.append(_extname('SPECTRUM'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(_detnam(detnam))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(date_obs))
    header.append(_dateend(date_end))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(tstart))
    header.append(_tstop(tstop))
    header.append(_trigtime(trigtime))
    header.append(('datatype', datatype, 'GBM datatype used for this file'))
    header.append(_object(object))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(ra_obj))
    header.append(_decobj(dec_obj))
    header.append(_errrad(err_rad))
    header.append(_filter())
    header.append(
        ('areascal', 1., 'No special scaling of effective area by channel'))
    header.append(('backfile', backfile,
                   'Name of corresponding background file (if any)'))
    header.append(('backscal', 1., 'background file scaling factor'))
    header.append(
        ('corrfile', 'none', 'Name of corresponding correction file (if any)'))
    header.append(('corrscal', 1., 'Correction scaling file'))
    header.append(
        ('respfile', rspfile, 'Name of corresponding RMF file (if any)'))
    header.append(
        ('ancrfile', 'none', 'Name of corresponding ARF file (if any)'))
    header.append(('sys_err', 0., 'No systematic errors'))
    header.append(('poisserr', poisserr, 'Assume Poisson Errors'))
    header.append(('grouping', 0, 'No special grouping has been applied'))
    header.append(_hduclass())
    header.append(
        ('hduclas1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'))
    header.append(
        ('hduclas2', 'TOTAL', 'Indicates gross data (source + background)'))
    header.append(('hduclas3', 'COUNT', 'Indicates data stored as counts'))
    header.append(('hduclas4', 'TYPEI', 'Indicates PHA Type I file format'))
    header.append(_hduvers())
    header.append(_chantype())
    header.append(_detchans(detchans))
    header.append(('exposure', exposure, 'Accumulation time - deadtime'))
    header.append(_extver())

    return header


def healpix_primary(tcat=None, trigtime=0.0):
    """Write the primary header of a HEALPix FITS file
    
    Parameters:
    -----------
    tcat: Tcat, optional
        The tcat object.  If set, then it will copy the relevant info from 
        the tcat into the primary header of the FITS file
    
    Returns:
    --------
    header: astropy.io.fits.header
        The header
    """

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(_creator(gbm.__version__))
    header.append(('filetype', 'IMAGE', 'Name for this type of FITS file'))
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(_date(current_time))
    header.append(_dateobs(''))
    header.append(_dateend(''))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(_tstart(0.0))
    header.append(_tstop(0.0))
    header.append(('filename', '', 'Name of this file'))
    header.append(_trigtime(trigtime))
    header.append(_object(''))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(_raobj(0.0))
    header.append(_decobj(0.0))
    header.append(_errrad(0.0))
    header.append(('theta', 0.0, '[deg] Angle from spacecraft zenith'))
    header.append(
        ('phi', 0.0, '[deg] Angle from spacecraft +X axis toward +Y'))
    header.append(('loc_src', 'Fermi, GBM',
                   'Mission/Instrument providing the localization'))
    header.append(('class', 'GRB', 'Classification of trigger'))
    header.append(('obj_clas', 'GRB', 'Classification of trigger'))
    header.append(
        ('geo_long', 0.0, '[deg] Spacecraft geographical east longitude'))
    header.append(
        ('geo_lat', 0.0, '[deg] Spacecraft geographical north latitude'))
    header.append(('ra_scx', 0.0, '[deg] Pointing of spacecraft x-axis: RA'))
    header.append(('dec_scx', 0.0, '[deg] Pointing of spacecraft x-axis: Dec'))
    header.append(('ra_scz', 0.0, '[deg] Pointing of spacecraft z-axis: RA'))
    header.append(('dec_scz', 0.0, '[deg] Pointing of spacecraft z-axis: Dec'))
    header.append(('loc_ver', '', 'Version string of localizing software'))
    header.append(
        ('loc_enrg', '(50, 300)', 'Energy range used for localization'))
    header.append(('lmethod', 'Interactive', 'Method of localization'))

    if tcat is not None:

        def insert_key(key):
            if key in tcat.headers['PRIMARY']:
                header[key] = tcat.headers['PRIMARY'][key]

        # copy tcat values to primary header
        insert_key('DATE-OBS')
        insert_key('DATE-END')
        insert_key('TSTART')
        insert_key('TSTOP')
        insert_key('TRIGTIME')
        insert_key('OBJECT')
        insert_key('RA_OBJ')
        insert_key('DEC_OBJ')
        insert_key('ERR_RAD')
        insert_key('THETA')
        insert_key('PHI')
        insert_key('GEO_LONG')
        insert_key('GEO_LAT')
        insert_key('RA_SCX')
        insert_key('DEC_SCX')
        insert_key('RA_SCZ')
        insert_key('DEC_SCZ')
        insert_key('LOC_VER')

    return header


def healpix_image(nside=128, object=None, extra_keys=None):
    """Write the image extension header of a HEALPix FITS file
    
    Parameters:
    -----------
    nside: int, optional
        The nside of the HEALPix map
    extra_keys: list, optional
        An additional keys to be added to the header
    
    Returns:
    --------
    header: astropy.io.fits.header
        The header
    """
    from healpy import nside2npix
    header = fits.Header()
    header.append(
        ('TTYPE1', 'PROBABILITY', 'Differential probability per pixel'))
    header.append(('TTYPE2', 'SIGNIFICANCE', 'Integrated probability'))
    header.append(('PIXTYPE', 'HEALPIX', 'HEALPIX pixelisation'))
    header.append(
        ('ORDERING', 'NESTED', 'Pixel ordering scheme, either RING or NESTED'))
    header.append(
        ('COORDSYS', 'C', 'Ecliptic, Galactic or Celestial (equatorial)'))
    header.append(
        ('EXTNAME', 'HEALPIX', 'name of this binary table extension'))
    header.append(('NSIDE', nside, 'Resolution parameter of HEALPIX'))
    header.append(('FIRSTPIX', 0, 'First pixel # (0 based)'))
    header.append(('LASTPIX', nside2npix(nside), 'Last pixel # (0 based)'))
    header.append(('INDXSCHM', 'IMPLICIT', 'Indexing: IMPLICIT or EXPLICIT'))
    header.append(
        ('OBJECT', object, 'Sky coverage, either FULLSKY or PARTIAL'))

    if extra_keys is not None:
        for key in extra_keys:
            header.append(key)

    return header
