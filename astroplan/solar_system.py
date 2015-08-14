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
from astropy.utils.data import download_file

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

def _get_sky_coord(time, location, jplephem_de430_index):
    """
    Calculate position of moon relative to Earth by subtracting the
    vector pointing from the Earth-moon barycenter to the object
    indicated by ``jplephem_de430_index`` by
    the vector from the Earth-moon barycenter to the Earth
    """
    from jplephem.spk import SPK
    kernel = SPK.open(get_spk_file())
    # if jplephem_de430_index == 2:
    #     cartesian_position = (kernel[0,2].compute(time.jd) +
    #                           kernel[2,299].compute(time.jd) -
    #                           kernel[0,3].compute(time.jd) -
    #                           kernel[3,399].compute(time.jd))
    # else:
    #     cartesian_position = (kernel[0,jplephem_de430_index].compute(time.jd) -
    #                           kernel[0,3].compute(time.jd) -
    #                           kernel[3,399].compute(time.jd))
    # if jplephem_de430_index < 3:
    #     to_planet_ind = int(str(jplephem_de430_index)+'99')
    #     cartesian_position = (kernel[0,jplephem_de430_index].compute(time.jd) +
    #                           kernel[jplephem_de430_index,to_planet_ind].compute(time.jd) -
    #                           kernel[0,3].compute(time.jd) -
    #                           kernel[3,399].compute(time.jd))
    # else:
    cartesian_position = (kernel[0,jplephem_de430_index].compute(time.jd) -
                          kernel[0,3].compute(time.jd) -
                          kernel[3,399].compute(time.jd))
    x, y, z = cartesian_position*u.km

    # Convert to GCRS coordinates
    cartrep = CartesianRepresentation(x=x, y=y, z=z)
    return SkyCoord(cartrep, frame=GCRS(obstime=time,
                                        obsgeoloc=location))

def mercury(time, location):
    return _get_sky_coord(time, location, 1)

def venus(time, location):
    return _get_sky_coord(time, location, 2)

def mars(time, location):
    return _get_sky_coord(time, location, 4)

def jupiter(time, location):
    return _get_sky_coord(time, location, 5)

def saturn(time, location):
    return _get_sky_coord(time, location, 6)

def uranus(time, location):
    return _get_sky_coord(time, location, 7)

def neptune(time, location):
    return _get_sky_coord(time, location, 8)

def pluto(time, location):
    return _get_sky_coord(time, location, 9)
