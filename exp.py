from astroplan import(Observer, ExoplanetHostTarget,
                      ExoplanetInTransitConstraint, AtNightConstraint,
                      AltitudeConstraint, observability_table)
import astropy.units as u
from astropy.time import Time

apo = Observer.at_site('apo')
time_range = Time(['2015-06-01 00:00:00', '2015-09-01 00:00:00'])

planet_list = ['HD 209458 b', 'HD 189733 b', 'GJ 1214 b']
targets = [ExoplanetHostTarget.from_name(planet) for planet in planet_list]

constraints = [ExoplanetInTransitConstraint(buffer_time=1*u.hour),
               AltitudeConstraint(min=70*u.deg),
               AtNightConstraint.twilight_civil()]

table = observability_table(constraints, apo, targets, time_range=time_range)
print(table)
