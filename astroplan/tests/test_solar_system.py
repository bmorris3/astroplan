from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astroplan.solar_system import (mercury, venus, mars, jupiter, saturn,
                                    uranus, neptune, pluto)
from astropy.time import Time
from astroplan import get_site
import astropy.units as u
from numpy.testing import assert_allclose
import ephem

def print_pyephem_planet_coords():
    time = Time("1980-08-13 00:00:00")
    apo = get_site("APO")
    obs = ephem.Observer()
    obs.lat = apo.latitude.radian
    obs.lon = apo.longitude.radian
    obs.elevation = apo.height.value
    obs.date = time.datetime

    pe_merc = ephem.Mercury()
    pe_merc.compute(obs)

    pe_v = ephem.Venus()
    pe_v.compute(obs)

    pe_m = ephem.Mars()
    pe_m.compute(obs)


    pe_j = ephem.Jupiter()
    pe_j.compute(obs)

    pe_s = ephem.Saturn()
    pe_s.compute(obs)


    pe_u = ephem.Uranus()
    pe_u.compute(obs)

    pe_n = ephem.Neptune()
    pe_n.compute(obs)

    pe_p = ephem.Pluto()
    pe_p.compute(obs)

    print([[float(pe_merc.a_ra), float(pe_merc.a_dec)],
            [float(pe_v.a_ra), float(pe_v.a_dec)],
            [float(pe_m.a_ra), float(pe_m.a_dec)],
            [float(pe_j.a_ra), float(pe_j.a_dec)],
            [float(pe_s.a_ra), float(pe_s.a_dec)],
            [float(pe_u.a_ra), float(pe_u.a_dec)],
            [float(pe_n.a_ra), float(pe_n.a_dec)],
            [float(pe_p.a_ra), float(pe_p.a_dec)]])

def test_planet_coords():
    separation_tolerance = 1*u.degree
    time = Time("1980-08-13 00:00:00")
    apo = get_site("APO")
    ap_merc = mercury(time, apo)
    ap_v = venus(time, apo)
    ap_m = mars(time, apo)
    ap_j = jupiter(time, apo)
    ap_s = saturn(time, apo)
    ap_u = uranus(time, apo)
    ap_n = neptune(time, apo)
    ap_p = pluto(time, apo)

    get_ra_dec_list = lambda fixedtgt: [fixedtgt.ra.radian,
                                        fixedtgt.dec.radian]

    astroplan_coords_radians = list(map(get_ra_dec_list, [ap_merc, ap_v, ap_m,
                                                          ap_j, ap_s, ap_u,
                                                          ap_n, ap_p]))
    # Calculate these coords with print_pyephem_planet_coords()
    pyephem_coords_radians = [[2.2621934197955693, 0.33759276848151226],
                              [1.672337996581023, 0.338442409078051],
                              [3.4620617385136305, -0.1382275166770391],
                              [2.896853311868943, 0.12429850509762044],
                              [3.0853219164084273, 0.063953773789147],
                              [4.005575325256486, -0.3138789902104957],
                              [4.52950822189444, -0.3789041061785018],
                              [3.571485523833843, 0.13743336907204634]]

    assert_allclose(astroplan_coords_radians, pyephem_coords_radians,
                    atol=separation_tolerance.to(u.rad).value)