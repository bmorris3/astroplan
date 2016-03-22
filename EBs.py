from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astroplan import (Observer, FixedTarget, SecondaryBody,
                       PrimaryEclipseConstraint, SecondaryEclipseConstraint,
                       is_observable)

from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.time import Time
import numpy as np

apo = Observer.at_site("APO")

# Get catalog of binaries from the Slawson et al. 2011 via Vizier
Vizier.ROW_LIMIT = -1
binary_table = Vizier.get_catalogs('J/AJ/142/160/table3')[0]

def get_target_from_table(kic_number):
    """
    Get an `~astroplan.FixedTarget` with associated `~astroplan.Secondarybody`
    for KIC ``kic_number``.

    Parameters
    ----------
    kic_number : int
        KIC number of eclipsing binary

    Returns
    -------
    primary : `~astroplan.FixedTarget`
        A `~astroplan.Fixedtarget` object representing the primary with an
        associated `~astroplan.SecondaryBody` representing the secondary.
    """
    if not kic_number in binary_table['KIC'].data:
        raise ValueError('KIC {0} is not in the catalog'.format(kic_number))

    idx = np.argwhere(binary_table['KIC'].data == kic_number)[0][0]
    period = binary_table['Per'].data[idx]*u.day
    sum_scaled_radii = binary_table['R1_R2'].data[idx]
    t0 = Time(2400000+binary_table['BJD'][idx], format='jd')
    ra = binary_table['_RAJ2000'][idx]*u.deg
    dec = binary_table['_DEJ2000'][idx]*u.deg

    # This is a rough approximation that shouldn't be applied here, but it
    # roughly works for our three targets (Winn 2011, Eqn 19):
    eclipse_duration = 0.8*period*sum_scaled_radii/np.pi

    # B component with associated orbital parameters
    secondary = SecondaryBody(name="KIC {0:d} B".format(kic_number),
                              period=period,
                              eclipse_duration=eclipse_duration,
                              time_inferior_conjunction=t0)

    primary = FixedTarget(coord=SkyCoord(ra=ra, dec=dec),
                          name="KIC {0:d} A".format(kic_number),
                          secondaries=secondary)
    return primary

kic_numbers = [9652632, 8719897, 12418816]
targets = [get_target_from_table(kic_number) for kic_number in kic_numbers]
time_range = Time(['2015-10-01', '2015-11-15'])
primary_eclipses = is_observable(PrimaryEclipseConstraint(), apo, targets,
                          time_range=time_range)
secondary_eclipses = is_observable(SecondaryEclipseConstraint(), apo, targets,
                            time_range=time_range)
print(np.any(primary_eclipses))
print(np.any(secondary_eclipses))
