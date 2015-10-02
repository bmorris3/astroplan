# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.tests.helper import assert_quantity_allclose

# Package
from ..target import FixedTarget, EclipsingFixedTarget


def test_FixedTarget_from_name():
    """
    Check that resolving target names with the `SkyCoord.from_name` constructor
    to produce a `FixedTarget` accurately resolves the coordinates of Polaris.
    """

    # Resolve coordinates with SkyCoord.from_name classmethod
    polaris_from_name = FixedTarget.from_name('Polaris')
    polaris_from_name = FixedTarget.from_name('Polaris', name='Target 1')
    # Coordinates grabbed from SIMBAD
    polaris_from_SIMBAD = SkyCoord('02h31m49.09456s', '+89d15m50.7923s')

    # Make sure separation is small
    assert polaris_from_name.coord.separation(polaris_from_SIMBAD) < 1*u.arcsec


def test_FixedTarget_ra_dec():
    """
    Confirm that FixedTarget.ra and FixedTarget.dec are the same as the
    right ascension and declination stored in the FixedTarget.coord variable -
    which is a SkyCoord
    """

    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')
    vega = FixedTarget(vega_coords, name='Vega')
    assert vega.coord == vega_coords, 'Store coordinates directly'
    assert vega.coord.ra == vega_coords.ra == vega.ra, ('Retrieve RA from '
                                                        'SkyCoord')
    assert vega.coord.dec == vega_coords.dec == vega.dec, ('Retrieve Dec from '
                                                           'SkyCoord')

def test_eclipsing():
    # Source: http://exoplanets.org/detail/HD_189733_b
    coord = SkyCoord(ra=300.18213945*u.deg, dec=22.71085126*u.deg, frame='icrs')
    epoch = Time(2454279.436714, format='jd')
    period = 2.21857567*u.day
    duration = 0.0760*u.day

    time = Time('2015-10-02 05:48')
    hd189 = EclipsingFixedTarget(coord, name='HD 189733 b', epoch=epoch,
                                 period=period, duration=duration)
    next_eclipse = hd189.eclipse_time(Time.now(), which='next')
    previous_eclipse = hd189.eclipse_time(Time.now(), which='previous')
    nearest_eclipse = hd189.eclipse_time(Time.now(), which='nearest')

    # Check answers versus the Czech Astronomical Society's Exoplanet Transit
    # Database webpage results for T0(HJD)=2453988.80336, Per=2.2185733 d:
    ETD_previous_eclipse = Time(2457296.696, format='jd')
    ETD_next_eclipse = Time(2457298.915, format='jd')
    ETD_nearest_eclipse = (ETD_previous_eclipse if
                           abs(ETD_previous_eclipse - time) <
                           abs(ETD_next_eclipse - time) else ETD_next_eclipse)

    tolerance = 15*u.min
    assert abs(next_eclipse - ETD_next_eclipse) < tolerance
    assert abs(previous_eclipse - ETD_previous_eclipse) < tolerance
    assert abs(nearest_eclipse - ETD_nearest_eclipse) < tolerance
