from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astroplan.solar_system import (mercury, venus, mars, jupiter, saturn,
                                    uranus, neptune, pluto)
from astropy.time import Time
from astroplan import get_site
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

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

import ephem
obs = ephem.Observer()
obs.lat = apo.latitude.radian
obs.lon = apo.longitude.radian

pe_merc = ephem.Mercury()
pe_merc.compute(obs)
assert_quantity_allclose([ap_merc.ra, ap_merc.dec],
                         [float(pe_merc.ra)*u.rad, float(pe_merc.dec)*u.rad],
                         atol=separation_tolerance)

pe_v = ephem.Venus()
pe_v.compute(obs)
assert_quantity_allclose([ap_v.ra, ap_v.dec],
                         [float(pe_v.ra)*u.rad, float(pe_v.dec)*u.rad],
                         atol=separation_tolerance)

pe_m = ephem.Mars()
pe_m.compute(obs)
assert_quantity_allclose([ap_m.ra, ap_m.dec],
                         [float(pe_m.ra)*u.rad, float(pe_m.dec)*u.rad],
                         atol=separation_tolerance)

pe_j = ephem.Jupiter()
pe_j.compute(obs)
assert_quantity_allclose([ap_j.ra, ap_j.dec],
                         [float(pe_j.ra)*u.rad, float(pe_j.dec)*u.rad],
                         atol=separation_tolerance)

pe_s = ephem.Saturn()
pe_s.compute(obs)
assert_quantity_allclose([ap_s.ra, ap_s.dec],
                         [float(pe_s.ra)*u.rad, float(pe_s.dec)*u.rad],
                         atol=separation_tolerance)

pe_u = ephem.Uranus()
pe_u.compute(obs)
assert_quantity_allclose([ap_u.ra, ap_u.dec],
                         [float(pe_u.ra)*u.rad, float(pe_u.dec)*u.rad],
                         atol=separation_tolerance)

pe_n = ephem.Neptune()
pe_n.compute(obs)
assert_quantity_allclose([ap_n.ra, ap_n.dec],
                         [float(pe_n.ra)*u.rad, float(pe_n.dec)*u.rad],
                         atol=separation_tolerance)

pe_p = ephem.Pluto()
pe_p.compute(obs)
assert_quantity_allclose([ap_p.ra, ap_p.dec],
                         [float(pe_p.ra)*u.rad, float(pe_p.dec)*u.rad],
                         atol=separation_tolerance)