from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

from astroplan import (FixedTarget, SecondaryBody, is_observable, Observer,
                       PrimaryEclipseConstraint)

coord = SkyCoord.from_name('HD 189733')
# Source: http://exoplanets.org/detail/HD_189733_b
t0 = Time(2454279.436714, format='jd')
period = 2.21857567*u.day
eclipse_duration = 0.0760*u.day

hd189 = FixedTarget(coord, name='HD 189733')
hd189b = SecondaryBody(name='HD 189733 b', time_inferior_conjunction=t0,
                       period=period, eclipse_duration=eclipse_duration)

next_eclipse = hd189b.inferior_conjunction_time(Time.now(), which='next')
print(next_eclipse.iso)

obs = Observer.at_site("APO")
targets = [hd189b]
time_range1 = Time(['2015-09-30 00:00', '2015-10-01 06:00'])
time_range2 = Time(['2015-09-30 00:00', '2015-10-01 02:00'])
constraints = PrimaryEclipseConstraint()

is_in_transit_in_range = is_observable(constraints, obs, targets, time_range=time_range1)
not_in_transit =is_observable(constraints, obs, targets, time_range=time_range2)
assert np.any(is_in_transit_in_range)
assert not np.any(not_in_transit)



# from astroplan import SecondaryBody
# # Source: http://exoplanets.org/detail/HD_189733_b
# t0 = Time(2454279.436714, format='jd')
# period = 2.21857567*u.day
# eclipse_duration = 0.0760*u.day
# hd189 = SecondaryBody(name='HD 189733 b', time_inferior_conjunction=t0,
#                       period=period, eclipse_duration=eclipse_duration)
# egr = hd189.egress_time(Time.now(), which='next')
# ing = hd189.ingress_time(Time.now(), which='next')
# mid = hd189.inferior_conjunction_time(Time.now(), which='next')
#
# print(hd189.phase([ing, mid, egr]))