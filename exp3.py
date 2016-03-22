from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_sun

from astroplan import (FixedTarget, Observer, time_grid_from_range,
                       is_always_observable, AirmassConstraint)

vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707*u.deg, dec=8.20163837*u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067*u.deg,
                     dec=89.26410897*u.deg), name="Polaris")


subaru = Observer.at_site("Subaru")
targets = [vega, rigel, polaris]
target = vega

time = Time('2001-02-03 04:05:06')
time_ranges = [Time([time, time+1*u.hour]) + offset
               for offset in np.arange(0, 400, 100)*u.day]

# From test_compare_airmass_constraint_and_observer
for time_range in time_ranges:
    max_airmass = 2
    # Check if each target meets airmass constraint in using Observer
    # always_from_observer = [all([subaru.altaz(time, target).secz < max_airmass
    #                              for time in time_grid_from_range(time_range)])
    #                         for target in targets]
    always_from_observer = all([subaru.altaz(time, target).secz < max_airmass
                             for time in time_grid_from_range(time_range)])

    # Check if each target meets altitude constraints using
    # is_always_observable and AirmassConstraint
    always_from_constraint = is_always_observable(AirmassConstraint(max_airmass),
                                                  subaru, target,
                                                  time_range=time_range)

    print(always_from_observer, always_from_constraint)



