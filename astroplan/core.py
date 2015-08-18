# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 Angle, Latitude, Longitude,
                                 UnitSphericalRepresentation, SphericalRepresentation)

import astropy.units as u
import datetime
from astropy.time import Time
import pytz
import numpy as np
################################################################################
# TODO: Temporary solution to IERS tables problems
from astropy.utils.data import download_file
from astropy.utils import iers

iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL,
                                                      cache=True))
################################################################################


from astropy.extern.six import string_types
from .exceptions import TargetNeverUpWarning, TargetAlwaysUpWarning
from .sites import get_site
from .moon import get_moon, moon_illumination, moon_phase_angle
import warnings

from abc import ABCMeta, abstractmethod

import numpy as np

__all__ = ["Observer", "Target", "FixedTarget", "NonFixedTarget",
           "Constraint", "TimeWindow", "AltitudeRange",
           "AboveAirmass", "MAGIC_TIME"]

#__doctest_requires__ = {'*': ['scipy.integrate']}

MAGIC_TIME = Time(-999, format='jd')

def _generate_24hr_grid(t0, start, end, N, for_deriv=False):
    """
    Generate a nearly linearly spaced grid of time durations.

    The midpoints of these grid points will span times from ``t0``+``start``
    to ``t0``+``end``, including the end points, which is useful when taking
    numerical derivatives.

    Parameters
    ----------
    t0 : `~astropy.time.Time`
        Time queried for, grid will be built from or up to this time.

    start : float
        Number of days before/after ``t0`` to start the grid.

    end : float
        Number of days before/after ``t0`` to end the grid.

    N : int
        Number of grid points to generate

    for_deriv : bool
        Generate time series for taking numerical derivative (modify
        bounds)?

    Returns
    -------
    `~astropy.time.Time`
    """

    if for_deriv:
        time_grid = np.concatenate([[start - 1/(N-1)],
                                    np.linspace(start, end, N)[1:-1],
                                    [end + 1/(N-1)]])*u.day
    else:
        time_grid = np.linspace(start, end, N)*u.day

    return t0 + time_grid

def _target_is_vector(target):
    if hasattr(target, '__iter__'):
        return True
    else:
        return False

def list_FixedTarget_to_SkyCoord(list_of_FixedTargets):
    """
    Convert a list of `~astroplan.core.FixedTarget` objects to a vector
    `~astropy.coordinates.SkyCoord` object.

    Parameters
    ----------
    list_of_FixedTargets : list
        `~astroplan.core.FixedTarget` objects

    Returns
    -------
    sc : `~astropy.coordinates.SkyCoord`
    """
    coord_list = [target.coord for target in list_of_FixedTargets]
    sc = SkyCoord(SkyCoord(coord_list).data.represent_as(
                  UnitSphericalRepresentation),
                  representation=UnitSphericalRepresentation)
    return sc

class Observer(object):
    """
    A container class for information about an observer's location and
    environment.

    TODO: write this docstring
    """
    @u.quantity_input(elevation=u.m)
    def __init__(self, location=None, timezone='UTC', name=None, latitude=None,
                 longitude=None, elevation=0*u.m, pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        longitude : float, str, `~astropy.units.Quantity` (optional)
            The longitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Longitude` object.

        latitude : float, str, `~astropy.units.Quantity` (optional)
            The latitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Latitude` object.

        elevation : `~astropy.units.Quantity` (optional), default = 0 meters
            The elevation of the observing location, with respect to sea
            level. Defaults to sea level.

        pressure : `~astropy.units.Quantity` (optional)
            The ambient pressure. Defaults to zero (i.e. no atmosphere).

        relative_humidity : float (optional)
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity` (optional)
            The ambient temperature.

        timezone : str or `datetime.tzinfo` (optional)
            The local timezone to assume. If a string, it will be passed through
            `pytz.timezone()` to produce the timezone object.

        description : str (optional)
            A short description of the telescope, observatory or observing
            location.
        """

        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and
                                 longitude is not None):
            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) an instance of '
                            'astropy.coordinates.EarthLocation or (2) '
                            'latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude.')

        # Accept various timezone inputs, default to UTC
        if isinstance(timezone, datetime.tzinfo):
            self.timezone = timezone
        elif isinstance(timezone, string_types):
            self.timezone = pytz.timezone(timezone)
        else:
            raise TypeError('timezone keyword should be a string, or an '
                            'instance of datetime.tzinfo')

    def __repr__(self):
        """
        String representation of the `~astroplan.Observer` object.

        Examples
        --------

        >>> from astroplan import Observer
        >>> keck = Observer.at_site("Keck", timezone="US/Hawaii")
        >>> print(keck)                                    # doctest: +FLOAT_CMP
        <Observer: name='Keck',
            location (lon, lat, el)=(-155.478333333 deg, 19.8283333333 deg, 4160.0 m),
            timezone=<DstTzInfo 'US/Hawaii' LMT-1 day, 13:29:00 STD>>
        """
        class_name = self.__class__.__name__
        attr_names = ['name', 'location', 'timezone', 'pressure', 'temperature',
                      'relative_humidity']
        attr_values = [getattr(self, attr) for attr in attr_names]
        attributes_strings = []
        for name, value in zip(attr_names, attr_values):
            if value is not None:
                # Format location for easy readability
                if name == 'location':
                    formatted_loc = ["{} {}".format(i.value, i.unit)
                                     for i in value.to_geodetic()]
                    attributes_strings.append(
                        "{} (lon, lat, el)=({})".format(name,
                                                        ", ".join(formatted_loc)))
                else:
                    if name != 'name':
                        value = repr(value)
                    else:
                        value = "'{}'".format(value)
                    attributes_strings.append("{}={}".format(name, value))
        return "<{}: {}>".format(class_name, ",\n    ".join(attributes_strings))

    @classmethod
    def at_site(cls, site_name, **kwargs):
        """
        Initialize an `~astroplan.core.Observer` object with a site name.

        Extra keyword arguments are passed to the `~astroplan.core.Observer`
        constructor (see `~astroplan.core.Observer` for available keyword
        arguments).

        Parameters
        ----------
        site_name : str
            Observatory name, must be resolvable with
            `~astroplan.sites.get_site`.

        Returns
        -------
        `~astroplan.core.Observer`
            Observer object.

        Examples
        --------
        Initialize an observer at Kitt Peak National Observatory:

        >>> from astroplan import Observer
        >>> import astropy.units as u
        >>> kpno_generic = Observer.at_site('kpno')
        >>> kpno_today = Observer.at_site('kpno', pressure=1*u.bar, temperature=0*u.deg_C)
        """
        name = kwargs.pop('name', site_name)
        if 'location' in kwargs:
            raise ValueError("Location kwarg should not be used if "
                             "initializing an Observer with Observer.at_site()")
        return cls(location=get_site(site_name), name=name, **kwargs)

    def astropy_time_to_datetime(self, astropy_time):
        """
        Convert the `~astropy.time.Time` object ``astropy_time`` to a
        localized `~datetime.datetime` object.

        Timezones localized with `~pytz`.

        Parameters
        ----------
        astropy_time : `~astropy.time.Time`
            Scalar or list-like.

        Returns
        -------
        `~datetime.datetime`
            Localized datetime, where the timezone of the datetime is
            set by the ``timezone`` keyword argument of the
            `~astroplan.Observer` constructor.
        """

        if not astropy_time.isscalar:
            return [self.astropy_time_to_datetime(t) for t in astropy_time]

        # Convert astropy.time.Time to a UTC localized datetime (aware)
        utc_datetime = pytz.utc.localize(astropy_time.utc.datetime)

        # Convert UTC to local timezone
        return self.timezone.normalize(utc_datetime)

    def datetime_to_astropy_time(self, date_time):
        """
        Convert the `~datetime.datetime` object ``date_time`` to a
        `~astropy.time.Time` object.

        Timezones localized with `~pytz`. If the ``date_time`` is naive, the
        implied timezone is the ``timezone`` structure of ``self``.

        Parameters
        ----------
        date_time : `~datetime.datetime` or list-like

        Returns
        -------
        `~astropy.time.Time`
            Astropy time object (no timezone information preserved).
        """

        if hasattr(date_time, '__iter__'):
            return Time([self.datetime_to_astropy_time(t) for t in date_time])

        # For timezone-naive datetimes, assign local timezone
        if date_time.tzinfo is None:
            date_time = self.timezone.localize(date_time)

        return Time(date_time, location=self.location)

    def _transform_target_list_to_altaz(self, times, targets):
        """
        Workaround for transforming a list of coordinates ``targets`` to
        altitudes and azimuths.

        Parameters
        ----------
        times : `~astropy.time.Time` or list of `~astropy.time.Time` objects
            Time of observation

        targets : `~astropy.coordinates.SkyCoord` or list of `~astropy.coordinates.SkyCoord` objects
            List of target coordinates

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        Returns
        -------
        altitudes : list
            List of altitudes for each target, at each time
        """
        if times.isscalar:
            times = Time([times])

        if not isinstance(targets, list) and targets.isscalar:
            targets = [targets]

        targets_is_unitsphericalrep = [x.data.__class__ is
                                       UnitSphericalRepresentation for x in targets]
        if all(targets_is_unitsphericalrep) or not any(targets_is_unitsphericalrep):
            repeated_times = np.tile(times, len(targets))
            ra_list = Longitude([x.icrs.ra for x in targets])
            dec_list = Latitude([x.icrs.dec for x in targets])
            repeated_ra = np.repeat(ra_list, len(times))
            repeated_dec = np.repeat(dec_list, len(times))
            inner_sc = SkyCoord(ra=repeated_ra, dec=repeated_dec)
            target_SkyCoord = SkyCoord(inner_sc.data.represent_as(UnitSphericalRepresentation),
                                       representation=UnitSphericalRepresentation)
            transformed_coord = target_SkyCoord.transform_to(AltAz(location=self.location,
                                                                   obstime=repeated_times))
        else:
            # TODO: This is super slow.
            repeated_times = np.tile(times, len(targets))
            repeated_targets = np.repeat(targets, len(times))
            target_SkyCoord = SkyCoord(SkyCoord(repeated_targets).data.represent_as(
                                       UnitSphericalRepresentation),
                                       representation=UnitSphericalRepresentation)

            transformed_coord = target_SkyCoord.transform_to(AltAz(location=self.location,
                                                                   obstime=repeated_times))
        return transformed_coord

    def altaz(self, time, target=None, obswl=None):
        """
        Get an `~astropy.coordinates.AltAz` frame or coordinate.

        If ``target`` is None, generates an altitude/azimuth frame. Otherwise,
        calculates the transformation to that frame for the requested ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            The time at which the observation is taking place. Will be used as
            the ``obstime`` attribute in the resulting frame or coordinate. This
            will be passed in as the first argument to the `~astropy.time.Time`
            initializer, so it can be anything that `~astropy.time.Time` will
            accept (including a `~astropy.time.Time` object)

        target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`, or list; defaults to `None` (optional)
            Celestial object(s) of interest. If ``target`` is `None`, return the
            `~astropy.coordinates.AltAz` frame without coordinates.

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation.

        Returns
        -------
        `~astropy.coordinates.AltAz`
            If ``target`` is `None`, returns `~astropy.coordinates.AltAz` frame.
            If ``target`` is not `None`, returns the ``target`` transformed to
            the `~astropy.coordinates.AltAz` frame.
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz_frame = AltAz(location=self.location, obstime=time,
                            pressure=self.pressure, obswl=obswl,
                            temperature=self.temperature,
                            relative_humidity=self.relative_humidity)
        if target is None:
            # Return just the frame
            return altaz_frame
        else:
            # If target is a list of targets:
            if _target_is_vector(target):
                get_coord = lambda x: x.coord if hasattr(x, 'coord') else x
                transformed_coords = self._transform_target_list_to_altaz(time,
                                          list(map(get_coord, target)))
                n_targets = len(target)
                new_shape = (n_targets, int(len(transformed_coords)/n_targets))

                for comp in transformed_coords.data.components:
                    getattr(transformed_coords.data, comp).resize(new_shape)
                return transformed_coords

            # If single target is a FixedTarget or a SkyCoord:
            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target
            return coordinate.transform_to(altaz_frame)

    def parallactic_angle(self, time, target):
        '''
        Calculate the parallactic angle.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Observation time.

        target : `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord` or list
            Target celestial object(s).

        Returns
        -------
        `~astropy.coordinates.Angle`
            Parallactic angle.

        Notes
        -----
        The parallactic angle is the angle between the great circle that
        intersects a celestial object and the zenith, and the object's hour
        circle [1]_.

        .. [1] https://en.wikipedia.org/wiki/Parallactic_angle

        '''
        if not isinstance(time, Time):
            time = Time(time)

        if _target_is_vector(target):
            get_coord = lambda x: x.coord if hasattr(x, 'coord') else x
            coordinate = SkyCoord(list(map(get_coord, target)))
        else:
            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target

        # Eqn (14.1) of Meeus' Astronomical Algorithms
        LST = time.sidereal_time('mean', longitude=self.location.longitude)
        H = (LST - coordinate.ra).radian
        q = np.arctan(np.sin(H) /
                      (np.tan(self.location.latitude.radian)*
                       np.cos(coordinate.dec.radian) -
                       np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad

        return Angle(q)

    # Sun-related methods.
    @u.quantity_input(horizon=u.deg)
    def _horiz_cross(self, t, alt, rise_set, horizon=0*u.degree):
        """
        Find time ``t`` when values in array ``a`` go from
        negative to positive or positive to negative (exclude endpoints)

        ``return_limits`` will return nearest times to zero-crossing.

        Parameters
        ----------
        t : `~astropy.time.Time`
            Grid of times
        alt : `~astropy.units.Quantity`
            Grid of altitudes
        rise_set : {"rising",  "setting"}
            Calculate either rising or setting across the horizon
        horizon : float
            Number of degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)
        Returns
        -------
        Returns the lower and upper limits on the time and altitudes
        of the horizon crossing.
        """
        alt = np.atleast_2d(Latitude(alt))
        n_targets = alt.shape[0]

        if rise_set == 'rising':
            # Find index where altitude goes from below to above horizon
            condition = (alt[:, :-1] < horizon) * (alt[:, 1:] > horizon)
        elif rise_set == 'setting':
            # Find index where altitude goes from above to below horizon
            condition = (alt[:, :-1] > horizon) * (alt[:, 1:] < horizon)

        target_inds, time_inds = np.nonzero(condition)

        if np.count_nonzero(condition) < n_targets:
            target_inds, _ = np.nonzero(condition)
            noncrossing_target_ind = np.setdiff1d(np.arange(n_targets),
                                                  target_inds,
                                                  assume_unique=True)#[0]

            warnmsg = ('Target(s) index {} does not cross horizon={} within '
                       '24 hours'.format(noncrossing_target_ind, horizon))

            if (alt[noncrossing_target_ind, :] > horizon).all():
                warnings.warn(warnmsg, TargetAlwaysUpWarning)
            else:
                warnings.warn(warnmsg, TargetNeverUpWarning)

            # Fill in missing time with MAGIC_TIME
            target_inds = np.insert(target_inds, noncrossing_target_ind,
                                    noncrossing_target_ind)
            time_inds = np.insert(time_inds.astype(float),
                                  noncrossing_target_ind,
                                  np.nan)
        elif np.count_nonzero(condition) > n_targets:
            old_target_inds = np.copy(target_inds)
            old_time_inds = np.copy(time_inds)

            time_inds = []
            target_inds = []
            for tgt, tm in zip(old_target_inds, old_time_inds):
                if tgt not in target_inds:
                    time_inds.append(tm)
                    target_inds.append(tgt)
            target_inds = np.array(target_inds)
            time_inds = np.array(time_inds)

        times = [t[i:i+2] if not np.isnan(i) else np.nan for i in time_inds]
        altitudes = [alt[i, j:j+2] if not np.isnan(j) else np.nan
                     for i, j in zip(target_inds, time_inds)]

        return times, altitudes

    @u.quantity_input(horizon=u.deg)
    def _two_point_interp(self, times, altitudes, horizon=0*u.deg):
        """
        Do linear interpolation between two ``altitudes`` at
        two ``times`` to determine the time where the altitude
        goes through zero.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Two times for linear interpolation between

        altitudes : array of `~astropy.units.Quantity`
            Two altitudes for linear interpolation between

        horizon : `~astropy.units.Quantity`
            Solve for the time when the altitude is equal to
            reference_alt.

        Returns
        -------
        t : `~astropy.time.Time`
            Time when target crosses the horizon

        """
        if not isinstance(times, Time):
            return MAGIC_TIME
        else:
            slope = (altitudes[1] - altitudes[0])/(times[1].jd - times[0].jd)
            return Time(times[1].jd - ((altitudes[1] - horizon)/slope).value,
                        format='jd')

    def _altitude_trig(self, LST, target):
        """
        Calculate the altitude of ``target`` at local sidereal times ``LST``.

        This method provides a factor of ~3 speed up over calling `altaz`, and
        inherently does *not* take the atmosphere into account.

        Parameters
        ----------
        LST : `~astropy.time.Time`
            Local sidereal times (array)

        target : {`~astropy.coordinates.SkyCoord`, `FixedTarget`} or similar
            Target celestial object's coordinates.

        Returns
        -------
        alt : `~astropy.unit.Quantity`
            Array of altitudes
        """
        alt = np.arcsin(np.sin(self.location.latitude.radian) *
                        np.sin(target.dec) +
                        np.cos(self.location.latitude.radian) *
                        np.cos(target.dec) *
                        np.cos(LST.radian - target.ra.radian))
        return alt

    def _calc_riseset(self, time, target, prev_next, rise_set, horizon, N=150):
        """
        Time at next rise/set of ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        rise_set : str - either 'rising' or 'setting'
            Compute prev/next rise or prev/next set

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        horizon : `~astropy.units.Quantity`
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of rise/set
        """

        if not isinstance(time, Time):
            time = Time(time)

        target_is_vector = _target_is_vector(target)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altaz = self.altaz(times, target)
        if target_is_vector:
            altitudes = [aa.alt for aa in altaz]
        else:
            altitudes = altaz.alt

        time_limits, altitude_limits = self._horiz_cross(times, altitudes, rise_set,
                                                    horizon)
        if not target_is_vector:
            return self._two_point_interp(time_limits[0], altitude_limits[0],
                                          horizon=horizon)
        else:
            return Time([self._two_point_interp(time_limit, altitude_limit,
                                                horizon=horizon)
                         for time_limit, altitude_limit in
                         zip(time_limits, altitude_limits)])

    def _calc_transit(self, time, target, prev_next, antitransit=False, N=150):
        """
        Time at next transit of the meridian of `target`.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        antitransit : bool
            Toggle compute antitransit (below horizon, equivalent to midnight
            for the Sun)

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of transit/antitransit
        """
        if not isinstance(time, Time):
            time = Time(time)

        target_is_vector = _target_is_vector(target)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N, for_deriv=True)
        else:
            times = _generate_24hr_grid(time, -1, 0, N, for_deriv=True)

        # The derivative of the altitude with respect to time is increasing
        # from negative to positive values at the anti-transit of the meridian
        if antitransit:
            rise_set = 'rising'
        else:
            rise_set = 'setting'

        altaz = self.altaz(times, target)
        if target_is_vector:
            d_altitudes = [each_alt.diff() for each_alt in altaz.alt]
        else:
            altitudes = altaz.alt
            d_altitudes = altitudes.diff()

        dt = Time((times.jd[1:] + times.jd[:-1])/2, format='jd')

        horizon = 0*u.degree # Find when derivative passes through zero
        time_limits, altitude_limits = self._horiz_cross(dt, d_altitudes,
                                                         rise_set, horizon)
        if not target_is_vector:
            return self._two_point_interp(time_limits[0], altitude_limits[0],
                                          horizon=horizon)
        else:
            return Time([self._two_point_interp(time_limit, altitude_limit,
                                                horizon=horizon)
                         for time_limit, altitude_limit in
                         zip(time_limits, altitude_limits)])

    def _determine_which_event(self, function, args_dict):
        """
        Run through the next/previous/nearest permutations of the solutions
        to `function(time, ...)`, and return the previous/next/nearest one
        specified by the args stored in args_dict.
        """
        time = args_dict.pop('time', None)
        target = args_dict.pop('target', None)
        which = args_dict.pop('which', None)
        horizon = args_dict.pop('horizon', None)
        rise_set = args_dict.pop('rise_set', None)
        antitransit = args_dict.pop('antitransit', None)

        # Assemble arguments for function, depending on the function.
        if function == self._calc_riseset:
            args = lambda w: (time, target, w, rise_set, horizon)
        elif function == self._calc_transit:
            args = lambda w: (time, target, w, antitransit)
        else:
            raise ValueError('Function {} not supported in '
                             '_determine_which_event.'.format(function))

        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_event = function(*args('next'))
            if which == 'next':
                return next_event

        if which == 'previous' or which == 'nearest':
            previous_event = function(*args('previous'))
            if which == 'previous':
                return previous_event

        if which == 'nearest':
            if _target_is_vector(target):
                return_times = []
                for next_e, prev_e in zip(next_event, previous_event):
                    if abs(time - prev_e) < abs(time - next_e):
                        return_times.append(prev_e)
                    else:
                        return_times.append(next_e)
                return Time(return_times)
            else:
                if abs(time - previous_event) < abs(time - next_event):
                    return previous_event
                else:
                    return next_event

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def target_rise_time(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate rise time.

        Compute time of the next/previous/nearest rise of the ``target``
        object, where "rise" is defined as the time when the ``target``
        transitions from altitudes below the ``horizon`` to above the
        ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Rise time of target
        """
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='rising',
                                                horizon=horizon))

    @u.quantity_input(horizon=u.deg)
    def target_set_time(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate set time.

        Compute time of the next/previous/nearest set of ``target``, where
        "set" is defined as when the ``target`` transitions from altitudes
        above ``horizon`` to below ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Set time of target.
        """
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='setting',
                                                horizon=horizon))

    def target_meridian_transit_time(self, time, target, which='nearest'):
        """
        Calculate time at the transit of the meridian.

        Compute time of the next/previous/nearest transit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Transit time of target
        """
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which,
                                                rise_set='setting'))

    def target_meridian_antitransit_time(self, time, target, which='nearest'):
        """
        Calculate time at the antitransit of the meridian.

        Compute time of the next/previous/nearest antitransit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Antitransit time of target
        """
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which, antitransit=True,
                                                rise_set='setting'))

    @u.quantity_input(horizon=u.deg)
    def sun_rise_time(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunrise.

        Compute time of the next/previous/nearest sunrise, where
        sunrise is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate.

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunrise
        """
        return self.target_rise_time(time, get_sun(time), which, horizon)

    @u.quantity_input(horizon=u.deg)
    def sun_set_time(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunset.

        Compute time of the next/previous/nearest sunset, where
        sunset is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.target_set_time(time, get_sun(time), which, horizon)

    def noon(self, time, which='nearest'):
        """
        Time at solar noon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar noon
        """
        return self.target_meridian_transit_time(time, get_sun(time), which)

    def midnight(self, time, which='nearest'):
        """
        Time at solar midnight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar midnight
        """
        return self.target_meridian_antitransit_time(time, get_sun(time), which)

    # Twilight convenience functions

    def twilight_evening_astronomical(self, time, which='nearest'):
        """
        Time at evening astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-18*u.degree)

    def twilight_evening_nautical(self, time, which='nearest'):

        """
        Time at evening nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-12*u.degree)

    def twilight_evening_civil(self, time, which='nearest'):
        """
        Time at evening civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-6*u.degree)

    def twilight_morning_astronomical(self, time, which='nearest'):
        """
        Time at morning astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_rise_time(time, which, horizon=-18*u.degree)

    def twilight_morning_nautical(self, time, which='nearest'):
        """
        Time at morning nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_rise_time(time, which, horizon=-12*u.degree)

    def twilight_morning_civil(self, time, which='nearest'):
        """
        Time at morning civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.sun_rise_time(time, which, horizon=-6*u.degree)

    # Moon-related methods.

    def moon_rise_time(self, time, **kwargs):
        """
        Returns the local moonrise time.

        The default moonrise returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_set_time(self, time, **kwargs):
        """
        Returns the local moonset time.

        The default moonset returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(self, time):
        """
        Calculate the illuminated fraction of the moon

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        moon : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the moon at time ``time``. If `None`, will calculate
            the position of the moon with `~astroplan.moon.get_moon`.

        sun : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the sun at time ``time``. If `None`, will calculate
            the position of the Sun with `~astropy.coordinates.get_sun`.

        Returns
        -------
        float
            Fraction of lunar surface illuminated
        """
        if not isinstance(time, Time):
            time = Time(time)

        return moon_illumination(time, self.location)

    def moon_phase(self, time=None, moon=None, sun=None):
        """
        Calculate lunar orbital phase.

        For example, phase=0 is "new", phase=1 is "full".

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        moon : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the moon at time ``time``. If `None`, will calculate
            the position of the moon with `~astroplan.moon.get_moon`.

        sun : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the sun at time ``time``. If `None`, will calculate
            the position of the Sun with `~astropy.coordinates.get_sun`.
        """
        if time is not None and not isinstance(time, Time):
            time = Time(time)

        return moon_phase_angle(time, self.location)

    def moon_altaz(self, time):
        """
        Returns the position of the moon in alt/az.

        TODO: Currently `moon_altaz` uses PyEphem to calculate the position
        of the moon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Returns
        -------
        altaz : `~astropy.coordinates.SkyCoord`
            Position of the moon transformed to altitude and azimuth
        """
        if not isinstance(time, Time):
            time = Time(time)

        try:
            import ephem
        except ImportError:
            raise ImportError("The moon_altaz function currently requires "
                              "PyEphem to compute the position of the moon.")

        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lat = self.location.latitude.to(u.degree).to_string(sep=':')
        obs.lon = self.location.longitude.to(u.degree).to_string(sep=':')
        obs.elevation = self.location.height.to(u.m).value
        if self.pressure is not None:
            obs.pressure = self.pressure.to(u.bar).value*1000.0

        if time.isscalar:
            obs.date = time.datetime
            moon.compute(obs)
            moon_alt = float(moon.alt)
            moon_az = float(moon.az)
        else:
            moon_alt = []
            moon_az = []
            for t in time:
                obs.date = t.datetime
                moon.compute(obs)
                moon_alt.append(float(moon.alt))
                moon_az.append(float(moon.az))
        return SkyCoord(alt=moon_alt*u.rad, az=moon_az*u.rad,
                        frame=self.altaz(time))

    @u.quantity_input(horizon=u.deg)
    def target_is_up(self, time, target, horizon=0*u.degree, return_altaz=False):
        """
        Is ``target`` above ``horizon`` at this ``time``?

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        return_altaz : bool (optional)
            Also return the '~astropy.coordinates.AltAz' coordinate.
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz = self.altaz(time, target)
        if _target_is_vector(target):
            observable = [alt > horizon for alt in altaz.alt]
        else:
            altitudes = altaz.alt
            observable = altitudes > horizon

        if not return_altaz:
            return observable
        else:
            return observable, altaz

    @u.quantity_input(horizon=u.deg)
    def is_night(self, time, horizon=0*u.deg, obswl=None):
        """
        Is the Sun below ``horizon`` at ``time``?

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating day/night (i.e.,
            -6 deg horizon = civil twilight, etc.)

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation

        Returns
        -------
        sun_below_horizon : bool
            `True` if sun is below ``horizon`` at ``time``, else `False`.
        """
        if not isinstance(time, Time):
            time = Time(time)

        solar_altitude = self.altaz(time, target=get_sun(time), obswl=obswl).alt
        return solar_altitude < horizon

class Target(object):
    """
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as ``FixedTarget`` or ``NonFixedTarget``.

    Would need to import six, abc to make this a metaclass?
    """
    __metaclass__ = ABCMeta

    def __init__(self, name=None, ra=None, dec=None, marker=None):
        """
        Defines a single observation target.

        Parameters
        ----------
        name : str, optional

        ra : WHAT TYPE IS ra ?

        dec : WHAT TYPE IS dec ?

        marker : str, optional
            User-defined markers to differentiate between different types
            of targets (e.g., guides, high-priority, etc.).
        """
        raise NotImplementedError()

    @property
    def ra(self):
        """
        Right ascension.
        """
        if isinstance(self, FixedTarget):
            return self.coord.ra
        if isinstance(self. NonFixedTarget):
            raise ValueError("NonFixedTarget objects require a time and/or "
                             "other specifications to calculate a position. "
                             "Did you mean to do NonFixedTarget.at(args).ra?")

    @property
    def dec(self):
        """
        Declination.
        """
        if isinstance(self, FixedTarget):
            return self.coord.dec
        if isinstance(self. NonFixedTarget):
            raise ValueError("NonFixedTarget objects require a time and/or "
                             "other specifications to calculate a position. "
                             "Did you mean to do NonFixedTarget.at(args).dec?")

class FixedTarget(Target):
    """
    An object that is "fixed" with respect to the celestial sphere.
    """
    def __init__(self, coord, name=None, **kwargs):
        """
        TODO: Docstring.
        """
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        """
        Initialize a `FixedTarget` by querying for a name, using the machinery
        in `~astropy.coordinates.SkyCoord.from_name`.
        """
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        if name is None:
            name = query_name
        return cls(SkyCoord.from_name(query_name), name=name, **kwargs)

    def __repr__(self):
        """
        String representation of `~astroplan.FixedTarget`.

        Examples
        --------
        Show string representation of a `~astroplan.FixedTarget` for Vega:

        >>> from astroplan import FixedTarget
        >>> from astroplan import FixedTarget
        >>> from astropy.coordinates import SkyCoord
        >>> vega_coord = SkyCoord(ra='279.23473479d', dec='38.78368896d')
        >>> vega = FixedTarget(coord=vega_coord, name="Vega")
        >>> print(vega)                             # doctest: +FLOAT_CMP
        <FixedTarget "Vega" at SkyCoord (ICRS): (ra, dec) in deg (279.23473479, 38.78368894)>
        """
        class_name = self.__class__.__name__
        fmt_coord = repr(self.coord).replace('\n   ', '')[1:-1]
        return '<{} "{}" at {}>'.format(class_name, self.name, fmt_coord)

class NonFixedTarget(Target):
    """
    An object that is not "fixed" with respect to the celestial sphere.

    Currently only some celestial objects are supported in this class.
    """
    def __init__(self, coord_function=None, name=None, constant_kwargs=None):
        """
        TODO: Docstring.
        """
        self.name = name.lower() if name is not None else name
        self.coord_function = coord_function
        self.constant_kwargs = constant_kwargs

    @classmethod
    def from_function(cls, coord_function, name=None):
        """
        Initialize a `~astropy.NonFixedTarget` by passing in a function that
        computes a `~astropy.coordinates.SkyCoord` for the target object.
        """
        return cls(coord_function=coord_function, name=name)

    def at(self, *args, **kwargs):
        """
        Get `~astropy.coordinates.SkyCoord` for the `~astroplan.NonFixedTarget`
        at a given time, location, etc.

        Parameters
        ----------
        All parameters passed to ``coord_function``.

        Returns
        -------
        A `~astropy.coordinates.SkyCoord` object as specified by the parameters.
        """
        if self.constant_kwargs is not None:
            for key in self.constant_kwargs:
                kwargs[key] = self.constant_kwargs[key]
        return self.coord_function(*args, **kwargs)

class Constraint(object):
    """
    An object containing observational constraints.

    A Constraints object is used in conjunction with a Target
    and an Observer object (via the apply_constraints method) to find out
    if a particular target is visible to the observer.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def apply_constraints(self, target, observer, constraint_list):
        """
        Returns information on a target's visibility.

        Finds out if a Target is observable by an Observer given a list
        of Constraint objects.  The list must contain at least one
        Constraint object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        constraint_list : WHAT TYPE IS constraint_list ? `numpy.array` ??
        """
        raise NotImplementedError


class TimeWindow(Constraint):
    """
    An object containing start and end times for an observation.
    """

    def __init__(self, start, end):
        """
        Initializes a TimeWindow object.

        Parameters
        ----------
        start : STRING OR astropy.time OBJECT ?

        end : STRING OR astropy.time OBJECT ?
        """
        raise NotImplementedError


class AltitudeRange(Constraint):
    """
    An object containing upper and lower altitude limits.
    """

    def __init__(self, low, high):
        """
        Initializes an AltitudeRange object.

        Parameters
        ----------
        low : `~astropy.units.Quantity`

        high : `~astropy.units.Quantity`
        """
        raise NotImplementedError


class AboveAirmass(Constraint):
    """
    An object containing an airmass lower limit.
    """

    def __init__(self, low):
        """
        Initializes an AboveAirmass object.

        Parameters
        ----------
        low : float
        """
        raise NotImplementedError
