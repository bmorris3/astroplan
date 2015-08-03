from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from ..exoplanets import Exoplanets
from ..core import Observer
from ..sites import get_site
from astropy.tests.helper import assert_quantity_allclose

def test_get_planet_index():
    db = Exoplanets()
    planet_name = 'HD 189733 b'
    idx = db.get_planet_index(planet_name)
    assert db.tbl['NAME'][idx] == planet_name

def test_get_ft_from_planet():
    db = Exoplanets()
    hd209 = "HD 209458 b"
    hd209_ft = db.get_skycoord(hd209)
    hd209_sc = SkyCoord(ra='22h03m10.77207s', dec='+18d53m03.5430s')

    gj1214 = "GJ 1214 b"
    gj1214_ft = db.get_skycoord(gj1214)
    gj1214_sc = SkyCoord(ra='17h15m18.94s', dec='+04d57m49.7s')

    assert_quantity_allclose([hd209_ft.ra, gj1214_ft.ra,
                              hd209_ft.dec, gj1214_ft.dec],
                             [hd209_sc.ra, gj1214_sc.ra,
                              hd209_sc.dec, gj1214_sc.dec],
                             atol=1*u.arcsec)

def test_database_parser():
    db = Exoplanets()
    planet = "HD 189733 b"
    idx = db.get_planet_index(planet)
    period_tbl = db.tbl['PER'][idx]
    period_exoplanetsorg = 2.21857567*u.day
    assert_quantity_allclose(period_tbl, period_exoplanetsorg,
                             atol=0.00001*u.day)

def test_get_visible_transits():
    # Verify with http://jefflcoughlin.com/transit.html or specifically
    # http://jefflcoughlin.com/cgi-bin/transit.cgi?Observatory=Kitt+Peak&Latitude=&Longitude=&JDstart=2456591.5&JDend=2456595.5&Sort=name&UTOpt=1&wavemu=1.0&albedo=0.0&fdist=0.6
    db = Exoplanets()
    start = Time(2456591.5, format='jd')
    end = Time(2456595.5, format='jd')

    kpno = Observer(location=get_site('kpno'))

    planet1 = 'CoRoT-1 b'
    visible_transits_corot1 = db.get_visible_transits(kpno, planet1, start, end)
    coughlin_corot1 = Time([2456591.889788, 2456594.907700], format='jd')
    assert all(abs(visible_transits_corot1 - coughlin_corot1) < 10*u.min)

    planet2 = 'WASP-50 b'
    visible_transits_wasp50 = db.get_visible_transits(kpno, planet2, start, end)
    coughlin_wasp50 = Time([2456592.857701, 2456594.812797], format='jd')
    assert all(abs(visible_transits_wasp50 - coughlin_wasp50) < 10*u.min)

