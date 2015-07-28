from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.constants import M_jup, M_sun, R_jup, R_sun
import numpy as np
from .core import FixedTarget

exoplanet_database_url = "http://www.exoplanets.org/csv-files/exoplanets.csv"
# TODO: Include a copy of the database in package? File size=~8MB
exoplanet_database_raw = download_file(exoplanet_database_url, cache=True)

exoplanet_table_units = dict(
    A=u.AU,
    AUPPER=u.AU,
    ALOWER=u.AU,
    UA=u.AU,
    AR=u.dimensionless_unscaled,
    ARUPPER=u.dimensionless_unscaled,
    ARLOWER=u.dimensionless_unscaled,
    UAR=u.dimensionless_unscaled,
    B=u.dimensionless_unscaled,
    BUPPER=u.dimensionless_unscaled,
    BLOWER=u.dimensionless_unscaled,
    UB=u.dimensionless_unscaled,
    BIGOM=u.deg,
    BIGOMUPPER=u.deg,
    BIGOMLOWER=u.deg,
    UBIGOM=u.deg,
    BMV=u.mag,
    CHI2=u.dimensionless_unscaled,
    DEC=u.deg,
    DENSITY=u.g/u.cm**3,
    DENSITYUPPER=u.g/u.cm**3,
    DENSITYLOWER=u.g/u.cm**3,
    UDENSITY=u.g/u.cm**3,
    DEPTH=u.dimensionless_unscaled,
    DEPTHUPPER=u.dimensionless_unscaled,
    DEPTHLOWER=u.dimensionless_unscaled,
    UDEPTH=u.dimensionless_unscaled,
    DIST=u.parsec,
    DISTUPPER=u.parsec,
    DISTLOWER=u.parsec,
    UDIST=u.parsec,
    DR=u.dimensionless_unscaled,
    DRUPPER=u.dimensionless_unscaled,
    DRLOWER=u.dimensionless_unscaled,
    UDR=u.dimensionless_unscaled,
    DVDT=u.m/u.s/u.day,
    DVDTUPPER=u.m/u.s/u.day,
    DVDTLOWER=u.m/u.s/u.day,
    UDVDT=u.m/u.s/u.day,
    ECC=u.dimensionless_unscaled,
    ECCUPPER=u.dimensionless_unscaled,
    ECCLOWER=u.dimensionless_unscaled,
    UECC=u.dimensionless_unscaled,
    FE=u.dimensionless_unscaled,
    FEUPPER=u.dimensionless_unscaled,
    FELOWER=u.dimensionless_unscaled,
    UFE=u.dimensionless_unscaled,
    GAMMA=u.km/u.s,
    GAMMAUPPER=u.km/u.s,
    GAMMALOWER=u.km/u.s,
    UGAMMA=u.km/u.s,
    GRAVITY=u.dex(u.cm/u.s**2),
    GRAVITYUPPER=u.dex(u.cm/u.s**2),
    GRAVITYLOWER=u.dex(u.cm/u.s**2),
    UGRAVITY=u.dex(u.cm/u.s**2),
    H=u.mag,
    I=u.deg,
    IUPPER=u.deg,
    ILOWER=u.deg,
    UI=u.deg,
    J=u.mag,
    K=u.m/u.s,
    KUPPER=u.m/u.s,
    KLOWER=u.m/u.s,
    UK=u.m/u.s,
    KS=u.mag,
    KP=u.mag,
    LAMBDA=u.deg,
    LAMBDAUPPER=u.deg,
    LAMBDALOWER=u.deg,
    ULAMBDA=u.deg,
    LOGG=u.deg,
    LOGGUPPER=u.deg,
    LOGGLOWER=u.deg,
    ULOGG=u.deg,
    MASS=M_jup,
    MASSUPPER=M_jup,
    MASSLOWER=M_jup,
    UMASS=M_jup,
    MSINI=M_jup,
    MSINIUPPER=M_jup,
    MSINILOWER=M_jup,
    UMSINI=M_jup,
    MSTAR=M_sun,
    MSTARUPPER=M_sun,
    MSTARLOWER=M_sun,
    UMSTAR=M_sun,
    OM=u.deg,
    OMUPPER=u.deg,
    OMLOWER=u.deg,
    UOM=u.deg,
    PAR=u.marcsec,
    PARUPPER=u.marcsec,
    PARLOWER=u.marcsec,
    UPAR=u.marcsec,
    PER=u.day,
    PERUPPER=u.day,
    PERLOWER=u.day,
    UPER=u.day,
    R=R_jup,
    RUPPER=R_jup,
    RLOWER=R_jup,
    UR=R_jup,
    RA=u.hourangle,
    RHK=u.dimensionless_unscaled,
    RHOSTAR=u.g/u.cm**3,
    RHOSTARUPPER=u.g/u.cm**3,
    RHOSTARLOWER=u.g/u.cm**3,
    URHOSTAR=u.g/u.cm**3,
    RMS=u.dimensionless_unscaled,
    RR=u.dimensionless_unscaled,
    RRUPPER=u.dimensionless_unscaled,
    RRLOWER=u.dimensionless_unscaled,
    URR=u.dimensionless_unscaled,
    RSTAR=R_sun,
    RSTARUPPER=R_sun,
    RSTARLOWER=R_sun,
    URSTAR=R_sun,
    SEDEPTHJ=u.dimensionless_unscaled,
    SEDEPTHJUPPER=u.dimensionless_unscaled,
    SEDEPTHJLOWER=u.dimensionless_unscaled,
    USEDEPTHJ=u.dimensionless_unscaled,
    SEDEPTHH=u.dimensionless_unscaled,
    SEDEPTHHUPPER=u.dimensionless_unscaled,
    SEDEPTHHLOWER=u.dimensionless_unscaled,
    USEDEPTHH=u.dimensionless_unscaled,
    SEDEPTHKS=u.dimensionless_unscaled,
    SEDEPTHKSUPPER=u.dimensionless_unscaled,
    SEDEPTHKSLOWER=u.dimensionless_unscaled,
    USEDEPTHKS=u.dimensionless_unscaled,
    SEDEPTHKP=u.dimensionless_unscaled,
    SEDEPTHKPUPPER=u.dimensionless_unscaled,
    SEDEPTHKPLOWER=u.dimensionless_unscaled,
    USEDEPTHKP=u.dimensionless_unscaled,
    SEDEPTH36=u.dimensionless_unscaled,
    SEDEPTH36UPPER=u.dimensionless_unscaled,
    SEDEPTH36LOWER=u.dimensionless_unscaled,
    USEDEPTH36=u.dimensionless_unscaled,
    SEDEPTH45=u.dimensionless_unscaled,
    SEDEPTH45UPPER=u.dimensionless_unscaled,
    SEDEPTH45LOWER=u.dimensionless_unscaled,
    USEDEPTH45=u.dimensionless_unscaled,
    SEDEPTH58=u.dimensionless_unscaled,
    SEDEPTH58UPPER=u.dimensionless_unscaled,
    SEDEPTH58LOWER=u.dimensionless_unscaled,
    USEDEPTH58=u.dimensionless_unscaled,
    SEDEPTH80=u.dimensionless_unscaled,
    SEDEPTH80UPPER=u.dimensionless_unscaled,
    SEDEPTH80LOWER=u.dimensionless_unscaled,
    USEDEPTH80=u.dimensionless_unscaled,
    SEP=u.AU,
    SEPUPPER=u.AU,
    SEPLOWER=u.AU,
    USEP=u.AU,
    SHK=u.dimensionless_unscaled,
    T0=u.day,
    T0UPPER=u.day,
    T0LOWER=u.day,
    UT0=u.day,
    T14=u.day,
    T14UPPER=u.day,
    T14LOWER=u.day,
    UT14=u.day,
    TEFF=u.Kelvin,
    TEFFUPPER=u.Kelvin,
    TEFFLOWER=u.Kelvin,
    UTEFF=u.Kelvin,
    TT=u.day,
    TTUPPER=u.day,
    TTLOWER=u.day,
    UTT=u.day,
    V=u.mag,
    VSINI=u.km/u.s,
    VSINIUPPER=u.km/u.s,
    VSINILOWER=u.km/u.s,
    UVSINI=u.km/u.s,
)

def parse_raw_database(raw):
    table = ascii.read(raw)
    for column in table.columns:
        if column in exoplanet_table_units:
            table[column].unit = exoplanet_table_units[column]
    return table

exoplanet_table = parse_raw_database(exoplanet_database_raw)

def get_planet_index(planet, tbl=exoplanet_table):
    return np.argwhere(planet == tbl['NAME'])[0][0]

def get_FixedTarget_from_exoplanet(planet, tbl=exoplanet_table):
    idx = get_planet_index(planet)
    sc = SkyCoord(ra=tbl['RA'].quantity[idx], dec=tbl['DEC'].quantity[idx])
    return FixedTarget(coord=sc, name=planet)

def get_transits(planet, time_start, time_end, tbl=exoplanet_table):
    """
    Get mid-transit times of ``planet`` between ``time_start`` and ``time_end``.

    Parameters
    ----------
    planet : str
        Planet name

    time_start : `~astropy.time.Time`
        Calculate transit times after ``time_start``

    time_end : `~astropy.time.Time`
        Calculate transit times before ``time_end``

    Returns
    -------
    times : `~astropy.time.Time`
        Transit times between ``time_start`` and ``time_end``
    """
    idx = get_planet_index(planet, tbl=tbl)
    period = tbl['PER'].quantity[idx]
    epoch = Time(tbl['TT'].quantity[idx], format='jd')
    start = np.ceil((time_start - epoch)/period)
    end = np.floor((time_end - epoch)/period)

    if start == end:
        return start*period + epoch
    elif start < end:
        return np.arange(start, end+1, dtype=int)*period + epoch
    else:
        return []

@u.quantity_input(horizon=u.deg)
def transit_visible(time, observatory, planet, horizon=0*u.degree):
    """
    Is transit of ``planet`` visible during night of ``time`` visible
    for an observer at ``observatory``?

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    observatory : `~astroplan.core.Observer`
        Observer

    planet : str
        Name of exoplanet
    """
    target = get_FixedTarget_from_exoplanet(planet)
    if observatory.can_see(time, target):
        rise_time = observatory.calc_rise(time, target, horizon=horizon,
                                          which='previous')
        set_time = observatory.calc_set(time, target, horizon=horizon,
                                        which='next')
    else:
        rise_time = observatory.calc_rise(time, target, which='next')
        set_time = observatory.calc_set(time, target, which='next')

    transits = get_transits(planet, rise_time, set_time)
    if hasattr(transits, 'isscalar') and transits.isscalar:
        return True
    elif hasattr(transits, 'isscalar') and not transits.isscalar:
        return len(transits) > 0
    else:
        return False
