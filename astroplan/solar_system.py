# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for interacting with jplephem to compute positions of solar system
objects.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (SkyCoord, GCRS, CartesianRepresentation)
import astropy.units as u
from astropy.utils.data import download_file
import numpy as np
from astropy.constants import c as speed_of_light
from jplephem.spk import SPK
from .core import NonFixedTarget

__all__ = ["mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune",
           "pluto"]

def _get_spk_file():
    """
    Get the Satellite Planet Kernel (SPK) file `de430.bsp` from NASA JPL.

    Download the file from the JPL webpage once and subsequently access a
    cached copy. This file is ~120 MB, and covers years ~1550-2650 CE [1]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    """
    de430_url = ('http://naif.jpl.nasa.gov/pub/naif/'
                 'generic_kernels/spk/planets/de430.bsp')
    return download_file(de430_url, cache=True, show_progress=True)

kernel = SPK.open(_get_spk_file())

def _get_cartesian_position(time, planet_index):
    """
    Calculate the position of planet ``planet_index`` in cartesian coordinates.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    cartesian_position : `~astropy.units.Quantity`
        GRCS position of the planet in cartesian coordinates
    """
    if planet_index < 3:
        barycenter_to_planet_ind = planet_index*100 + 99
        cartesian_position = (kernel[0,planet_index].compute(time.jd) +
                              kernel[planet_index,barycenter_to_planet_ind].compute(time.jd) -
                              kernel[0,3].compute(time.jd) -
                              kernel[3,399].compute(time.jd))
    else:
        cartesian_position = (kernel[0,planet_index].compute(time.jd) -
                              kernel[0,3].compute(time.jd) -
                              kernel[3,399].compute(time.jd))
    return cartesian_position*u.km

def _get_earth_distance(time, planet_index):
    """
    Calculate the distance between Earth and planet ``planet_index``.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and planet.
    """
    if planet_index < 3:
        barycenter_to_planet_ind = planet_index*100 + 99
        cartesian_position_planet = (kernel[0,planet_index].compute(time.jd) +
                                     kernel[planet_index,barycenter_to_planet_ind].compute(time.jd))
    else:
        cartesian_position_planet = kernel[0,planet_index].compute(time.jd)

    cartesian_position_earth = (kernel[0,3].compute(time.jd) +
                                kernel[3,399].compute(time.jd))

    # Quadrature sum of the difference of the position vectors
    earth_distance = np.sqrt(np.sum((np.array(cartesian_position_planet) -
                                     np.array(cartesian_position_earth))**2))*u.km

    return earth_distance

def _get_sky_coord(time, location, planet_index):
    """
    Create a `~astropy.coordinates.SkyCoord` object for planet ``planet_index``.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the planet
    """
    # Get distance of planet at `time`
    earth_distance = _get_earth_distance(time, planet_index)
    light_travel_time = earth_distance/speed_of_light
    emitted_time = time - light_travel_time

    # Calculate position given approximate light travel time
    x, y, z = _get_cartesian_position(emitted_time, planet_index)
    cartrep = CartesianRepresentation(x=x, y=y, z=z)
    return SkyCoord(cartrep, frame=GCRS(obstime=time,
                                        obsgeoloc=location))

def mercury(location=None):
    """
    Position of the planet Mercury.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation` (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Mercury",
                          constant_kwargs=dict(location=location,
                                               planet_index=1))

def venus(location=None):
    """
    Position of the planet Venus.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`  (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Venus",
                          constant_kwargs=dict(location=location,
                                               planet_index=2))

def mars(location=None):
    """
    Position of the planet Mars.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`  (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Mars",
                          constant_kwargs=dict(location=location,
                                               planet_index=4))

def jupiter(location=None):
    """
    Position of the planet Jupiter.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`  (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Jupiter",
                          constant_kwargs=dict(location=location,
                                               planet_index=5))

def saturn(location):
    """
    Position of the planet Saturn.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`  (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Saturn",
                          constant_kwargs=dict(location=location,
                                               planet_index=6))

def uranus(location=None):
    """
    Position of the planet Uranus.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Uranus",
                          constant_kwargs=dict(location=location,
                                               planet_index=7))

def neptune(location=None):
    """
    Position of the planet Neptune.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`  (optional)
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Neptune",
                          constant_kwargs=dict(location=location,
                                               planet_index=8))

def pluto(location=None):
    """
    Position of the planet Pluto.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return NonFixedTarget(coord_function=_get_sky_coord,
                          name="Pluto",
                          constant_kwargs=dict(location=location,
                                               planet_index=9))
