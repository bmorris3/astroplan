# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for interacting with jplephem to compute positions of solar system
objects.

DE430 key:
    File type DAF/SPK and format LTL-IEEE with 14 segments:
    Solar System Barycenter (0) -> Mercury Barycenter (1)
    Solar System Barycenter (0) -> Venus Barycenter (2)
    Solar System Barycenter (0) -> Earth Barycenter (3)
    Solar System Barycenter (0) -> Mars Barycenter (4)
    Solar System Barycenter (0) -> Jupiter Barycenter (5)
    Solar System Barycenter (0) -> Saturn Barycenter (6)
    Solar System Barycenter (0) -> Uranus Barycenter (7)
    Solar System Barycenter (0) -> Neptune Barycenter (8)
    Solar System Barycenter (0) -> Pluto Barycenter (9)
    Solar System Barycenter (0) -> Sun (10)
    Earth Barycenter (3) -> Moon (301)
    Earth Barycenter (3) -> Earth (399)
    Mercury Barycenter (1) -> Mercury (199)
    Venus Barycenter (2) -> Venus (299)
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

def get_spk_file():
    """
    Get the Satellite Planet Kernel (SPK) file `de430.bsp` from NASA JPL.

    Download the file from the JPL webpage once and subsequently access a
    cached copy. This file is ~120 MB, and covers years ~1550-2650 CE [1]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    """
    de430_url = ('http://naif.jpl.nasa.gov/pub/naif/'
                 'generic_kernels/spk/planets/de430.bsp')
    return download_file(de430_url, cache=True, show_progress=True)

kernel = SPK.open(get_spk_file())

def _get_cartesian_position(time, planet_index):
    """
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

def mercury(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 1),
                       name="Mercury")

def venus(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 2),
                       name="Venus")

def mars(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 4),
                       name="Mars")

def jupiter(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 5),
                       name="Jupiter")

def saturn(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 6),
                       name="Saturn")

def uranus(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 7),
                       name="Uranus")

def neptune(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 8),
                       name="Neptune")

def pluto(time, location):
    return FixedTarget(coord=_get_sky_coord(time, location, 9),
                       name="Pluto")
