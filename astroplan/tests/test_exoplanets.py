from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from ..exoplanets import (exoplanet_table, get_FixedTarget_from_exoplanet, \
                          get_planet_index, get_transits)
from astropy.tests.helper import assert_quantity_allclose

def test_get_planet_index():
    planet_name = 'HD 189733 b'
    idx = get_planet_index(planet_name)
    assert exoplanet_table['NAME'][idx] == planet_name

def test_get_ft_from_planet():
    hd209 = "HD 209458 b"
    hd209_ft = get_FixedTarget_from_exoplanet(hd209)
    hd209_sc = SkyCoord(ra='22h03m10.77207s', dec='+18d53m03.5430s')

    gj1214 = "GJ 1214 b"
    gj1214_ft = get_FixedTarget_from_exoplanet(gj1214)
    gj1214_sc = SkyCoord(ra='17h15m18.94s', dec='+04d57m49.7s')

    assert_quantity_allclose([hd209_ft.ra, gj1214_ft.ra,
                              hd209_ft.dec, gj1214_ft.dec],
                             [hd209_sc.ra, gj1214_sc.ra,
                              hd209_sc.dec, gj1214_sc.dec],
                             atol=1*u.arcsec)
def test_database_parser():
    planet = "HD 189733 b"
    idx = get_planet_index(planet)
    period_tbl = exoplanet_table['PER'].quantity[idx]
    period_exoplanetsorg = 2.21857567*u.day
    assert_quantity_allclose(period_tbl, period_exoplanetsorg,
                             atol=0.00001*u.day)

