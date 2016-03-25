# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for interacting with jplephem to compute positions of solar system
objects.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (SkyCoord, GCRS, CartesianRepresentation)
import astropy.units as u
from .core import FixedTarget
from astropy.utils.data import download_file
import numpy as np
from astropy.constants import c as speed_of_light
from jplephem.spk import SPK

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


def _get_absolute_planet_position(time, planet_index):
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

    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and planet.
    Notes
    -----

    """
    if planet_index < 3:
        barycenter_to_planet_ind = planet_index*100 + 99
        cartesian_position_planet = (kernel[0, planet_index].compute(time.jd) +
                                     kernel[planet_index,
                                            barycenter_to_planet_ind].compute(time.jd))
    else:
        cartesian_position_planet = kernel[0, planet_index].compute(time.jd)

    cartesian_position_earth = (kernel[0, 3].compute(time.jd) +
                                kernel[3, 399].compute(time.jd))

    # Quadrature sum of the difference of the position vectors
    earth_distance = np.sqrt(np.sum((np.array(cartesian_position_planet) -
                                     np.array(cartesian_position_earth))**2))*u.km

    earth_to_planet_vector = u.Quantity(kernel[0, planet_index].compute(time.jd) -
                                        cartesian_position_earth,
                                        unit=u.km)

    return earth_to_planet_vector, earth_distance


def _get_apparent_planet_position(time, planet_index):
    """
    Calculate the apparent position of planet ``planet_index`` in cartesian
    coordinates, given the approximate light travel time to the object.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Position of the planet defined as a vector from the Earth.

    Notes
    -----

    """
    # Get distance of planet at `time`
    earth_to_planet_vector, earth_distance = _get_absolute_planet_position(time, planet_index)

    # The apparent position depends on the time that the light was emitted from
    # the distant planet, so subtract off the light travel time
    light_travel_time = earth_distance/speed_of_light
    emitted_time = time - light_travel_time

    # Calculate position given approximate light travel time.
    # TODO: this should be solved iteratively to converge on precise positions
    earth_to_planet_vector, earth_distance = _get_absolute_planet_position(emitted_time, planet_index)
    x, y, z = earth_to_planet_vector
    return CartesianRepresentation(x=x, y=y, z=z)


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
    cartrep = _get_apparent_planet_position(time, planet_index)

    return SkyCoord(GCRS(cartrep, obstime=time,
                         obsgeoloc=u.Quantity(location.geocentric, copy=False)))


def mercury(time, location):
    """
    Position of the planet Mercury.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 1),
                       name="Mercury")


def venus(time, location):
    """
    Position of the planet Venus.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 2),
                       name="Venus")


def mars(time, location):
    """
    Position of the planet Mars.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 4),
                       name="Mars")


def jupiter(time, location):
    """
    Position of the planet Jupiter.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 5),
                       name="Jupiter")


def saturn(time, location):
    """
    Position of the planet Saturn.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 6),
                       name="Saturn")


def uranus(time, location):
    """
    Position of the planet Uranus.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 7),
                       name="Uranus")


def neptune(time, location):
    """
    Position of the planet Neptune.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 8),
                       name="Neptune")


def pluto(time, location):
    """
    Position of the planet Pluto.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth.

    Returns
    -------

    """
    return FixedTarget(coord=_get_sky_coord(time, location, 9),
                       name="Pluto")
