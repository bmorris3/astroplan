# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module relies on the exoplanet database at exoplanets.org for information
about exoplanets and their host stars, by Wright, J.T., Fakhouri, O., Marcy,
G.W., et al. 2011, PASP, 123, 412  [1]_.

.. [1] http://arxiv.org/abs/1409.7709
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.constants import M_jup, M_sun, R_jup, R_sun
import numpy as np
from astropy.table import QTable

__all__ = ["get_FixedTarget_from_exoplanet", "is_transit_visible",
           "get_visible_transits", "get_transits", "get_available_planets"]

exoplanet_database_url = "http://www.exoplanets.org/csv-files/exoplanets.csv"
# TODO: Include a copy of the database in package? File size=~8MB
exoplanet_database_raw = download_file(exoplanet_database_url, cache=True)

# Collect units associated with each numeric data column of the table
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
    DEC=u.degree,
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
    # TT=u.day,
    # TTUPPER=u.day,
    # TTLOWER=u.day,
    # UTT=u.day,
    V=u.mag,
    VSINI=u.km/u.s,
    VSINIUPPER=u.km/u.s,
    VSINILOWER=u.km/u.s,
    UVSINI=u.km/u.s,
)

def parse_raw_database(raw_CSV_file):
    """
    Read the exoplanets.org CSV database file [1]_ into a
    `~astropy.table.QTable`.

    Do some data whitening (i.e. convert times to JD), initialize astropy
    objects and apply units wherever possible.

    .. [1] http://www.exoplanets.org/csv-files/exoplanets.csv

    Parameters
    ----------
    raw_CSV_file : str
        Path to exoplanets.org CSV database

    Returns
    -------
    table : `~astropy.table.QTable`

    """
    table = ascii.read(raw_CSV_file, format='fast_csv')
    for column in table.columns:
        if column in exoplanet_table_units:
            table[column].unit = exoplanet_table_units[column]
    table = QTable(table)
    # Clean up mid-transit epoch:
    KOIs = np.array(["KOI" in name for name in table['NAME']])
    midtransit_epochs = np.array(table["TT"])
    midtransit_epochs[KOIs] += 2440000
    midtransit_epochs_astropy = []
    for i, t in enumerate(midtransit_epochs):
        if isinstance(t, float):
            midtransit_epochs_astropy.append(Time(t, format='jd'))
        else:
            midtransit_epochs_astropy.append(np.nan)

    table["TT_validated"] = midtransit_epochs_astropy

    # Using mixin column, this is an order of magnitude faster than
    # initializing a SkyCoord for each planet individually
    table["SkyCoord"] = SkyCoord(ra=table['RA'], dec=table['DEC'], frame='icrs')
    return table

def validate_transit_database(tbl):
    """
    Create a copy of the exoplanet database table ``tbl`` for transit
    computations.

    The new table will only have the planets that have well-defined periods and
    mid-transit epochs.

    Parameters
    ----------
    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.
    """
    valid_periods = tbl['PER'] != 0 *u.day
    valid_epochs = tbl['TT_validated'] != 0
    return tbl[valid_periods*valid_epochs]

complete_exoplanet_table = parse_raw_database(exoplanet_database_raw)
transiting_exoplanets = validate_transit_database(complete_exoplanet_table)

def get_planet_index(planet, tbl=complete_exoplanet_table):
    """
    Get the row index of ``planet`` in exoplanet database table ``tbl``.

    Parameters
    ----------
    planet : str
        Exoplanet name

    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.

    Returns
    -------
    index : int
        Row index of the exoplanet database table corresponding to ``planet``
    """
    index_array = np.argwhere(tbl['NAME'] == planet)
    if len(index_array) == 0:
        raise ValueError('Exoplanet "{}" not found.'.format(planet))
    return index_array[0][0]

def get_skycoord(planet, tbl=complete_exoplanet_table):
    """
    Get `~astropy.coordinates.SkyCoord` object for exoplanet ``planet``.

    Parameters
    ----------
    planet : str
        Name of exoplanet in exoplanets.org database [1]_.

    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.

    .. [1] http://exoplanets.org
    """
    idx = get_planet_index(planet, tbl=tbl)
    return tbl['SkyCoord'][idx]

def get_transits(planet, time_start, time_end, tbl=transiting_exoplanets):
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

    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.

    Returns
    -------
    times : `~astropy.time.Time`
        Transit times between ``time_start`` and ``time_end``
    """
    idx = get_planet_index(planet, tbl=tbl)
    period = tbl['PER'][idx]
    epoch = tbl['TT_validated'][idx]
    start = np.ceil((time_start - epoch)/period)
    end = np.floor((time_end - epoch)/period)

    if start == end:
        return start*period + epoch
    elif start < end:
        return Time(np.arange(start, end+1, dtype=int)*period + epoch)
    else:
        return None

@u.quantity_input(horizon=u.deg, solar_horizon=u.deg)
def is_transit_visible(observer, planet, mid_transit_time,
                       planet_horizon=0*u.degree, solar_horizon=-6*u.deg,
                       full_transit=True, tbl=transiting_exoplanets):
    """
    Is transit of ``planet`` at ``time`` visible at night for this ``observer``?

    "At night" is defined as times when the Sun is below ``solar_horizon``.

    Parameters
    ----------
    mid_transit_time : `~astropy.time.Time`
        Mid-transit time

    observer : `~astroplan.core.Observer`
        Observer

    planet : str
        Name of exoplanet

    planet_horizon : `~astropy.units.Quantity` (optional)
        Minimum altitude of the planet required to be "visible"

    solar_horizon : `~astropy.units.Quantity` (optional)
        Maximum altitude of the Sun defining which transits occur "at night"

    full_transit : bool
        Require visible from first to fourth contact (True) or only
        at mid-transit time (False)

    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.
    """
    if full_transit:
        idx = get_planet_index(planet, tbl=tbl)
        duration = tbl['T14'][idx]
        check_times = [mid_transit_time-duration, mid_transit_time+duration]
    else:
        check_times = [mid_transit_time]

    target = get_skycoord(planet)
    visible = [True if observer.can_see(t, target, horizon=planet_horizon) and
                       observer.is_night(t, horizon=solar_horizon)
               else False for t in check_times]
    return all(visible)

@u.quantity_input(planet_horizon=u.deg, solar_horizon=u.deg)
def get_visible_transits(observer, planet, start_time, end_time,
                         planet_horizon=0*u.degree, solar_horizon=-6*u.deg,
                         full_transit=True, tbl=transiting_exoplanets):
    """
    Calculate mid-transit times of visible transits within a time range.

    Parameters
    ----------
    observer : `~astroplan.core.Observer`
        Observer location and environment.

    planet : str
        Exoplanet name

    start_time : `~astropy.time.Time`
        Earliest time to consider

    end_time : `~astropy.time.Time`
        Latest time to consider

    planet_horizon : `~astropy.units.Quantity` (optional)
        Minimum altitude of the planet during the transit

    solar_horizon : `~astropy.units.Quantity` (optional)
        Maximum altitude of the Sun at the time of transit for the transit
        to be considered "visible"

    full_transit : bool (optional)
        If True, the transit will only be considered visible if the planet is
        visible from first to fourth contact (the full transit duration), else
        only checks the mid-transit time.

    tbl : `~astropy.table.QTable` (optional)
        Exoplanet database table.
    """
    all_transits = get_transits(planet, start_time, end_time, tbl=tbl)
    is_transit_visible_kwargs = dict(planet_horizon=planet_horizon,
                                     solar_horizon=solar_horizon,
                                     full_transit=full_transit, tbl=tbl)
    if all_transits is None:
        return None

    if all_transits.isscalar:
        if is_transit_visible(observer, planet, all_transits,
                              **is_transit_visible_kwargs):
            return all_transits
        else:
            return None
    else:
        print(all_transits, type(all_transits))
        visible_transits = [t for t in all_transits if (t is not None and
                            is_transit_visible(observer, planet, t,
                                               **is_transit_visible_kwargs))]
        if len(visible_transits) > 0:
            return Time(visible_transits)
        else:
            return None

def get_available_planets():
    """
    List of planets available in the local copy of the exoplanets.org database.

    Returns
    -------
    list
        Available planets
    """
    return list(complete_exoplanet_table['NAME'])
