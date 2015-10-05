# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
from abc import ABCMeta

# Third-party
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

__all__ = ["Target", "FixedTarget", "NonFixedTarget", "SecondaryBody"]

class Target(object):
    """
    Abstract base class for target objects.

    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `~astroplan.target.FixedTarget` or
    `~astroplan.target.NonFixedTarget`.
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
        raise NotImplementedError()

    @property
    def dec(self):
        """
        Declination.
        """
        if isinstance(self, FixedTarget):
            return self.coord.dec
        raise NotImplementedError()


class FixedTarget(Target):
    """
    Coordinates and metadata for an object that is "fixed" with respect to the
    celestial sphere.

    Examples
    --------
    Create a `~astroplan.FixedTarget` object for Sirius:

    >>> from astroplan import FixedTarget
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> sirius_coord = SkyCoord(ra=101.28715533*u.deg, dec=16.71611586*u.deg)
    >>> sirius = FixedTarget(coord=sirius_coord, name="Sirius")

    Create an equivalent `~astroplan.FixedTarget` object for Sirius by querying
    for the coordinates of Sirius by name:

    >>> from astroplan import FixedTarget
    >>> sirius = FixedTarget.from_name("Sirius")
    """
    def __init__(self, coord, name=None, secondaries=None, **kwargs):
        """
        Parameters
        ----------
        coord : `~astropy.coordinates.SkyCoord`
            Coordinate of the target

        name : str (optional)
            Name of the target, used for plotting and representing the target
            as a string

        secondaries : list or `None` (optional)
            List of `~astroplan.SecondaryBody` objects that orbit the target
        """
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord
        if not hasattr(secondaries, '__len__') and secondaries is not None:
            secondaries = [secondaries]
        self.secondaries = secondaries

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        """
        Initialize a `FixedTarget` by querying for a name, using the machinery
        in `~astropy.coordinates.SkyCoord.from_name`.

        Parameters
        ----------
        query_name : str
            Name of the target used to query for coordinates.

        name : string or `None`
            Name of the target to use within astroplan. If `None`, query_name
            is used as ``name``.

        Examples
        --------
        >>> from astroplan import FixedTarget
        >>> sirius = FixedTarget.from_name("Sirius")
        >>> sirius.coord                              # doctest: +FLOAT_CMP
        <SkyCoord (ICRS): (ra, dec) in deg
            (101.28715533, -16.71611586)>
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

    @classmethod
    def _from_name_mock(cls, query_name, name=None):
        """
        Mock method to replace `FixedTarget.from_name` in tests without
        internet connection.
        """
        stars = {
            "rigel": {"ra": 78.63446707*u.deg, "dec": -8.20163837*u.deg},
            "sirius": {"ra": 101.28715533*u.deg, "dec": -16.71611586*u.deg},
            "vega": {"ra": 279.23473479*u.deg, "dec": 38.78368896*u.deg},
            "aldebaran": {"ra": 68.98016279*u.deg, "dec": 16.50930235*u.deg},
            "polaris": {"ra": 37.95456067*u.deg, "dec": 89.26410897*u.deg}
        }

        if query_name.lower() in stars:
            return cls(coord=SkyCoord(**stars[query_name.lower()]),
                       name=query_name)
        else:
            raise ValueError("Target named {} not in mocked FixedTarget "
                             "method".format(query_name))

class SecondaryBody(object):
    """
    A body that orbits another body.
    """
    def __init__(self, name=None, period=None, time_inferior_conjunction=None,
                 eclipse_duration=None):
        """
        Parameters
        ----------
        name : str (optional)
            Name of the target, used for plotting and representing the target
            as a string

        time_inferior_conjunction : `~astropy.time.Time`
            Time at inferior conjunction (i.e. in eclipsing systems: the time
            of inferior conjunction is the time at the middle of primary
            eclipse)

        period : `~astropy.units.Quantity`
            Orbital period of the secondary

        eclipse_duration : `~astropy.units.Quantity` (optional)
            Duration of eclipse, if this secondary body eclipses the primary.
        """
        self.name = name
        self.period = period
        self.eclipse_duration = eclipse_duration
        self.time_inferior_conjunction = time_inferior_conjunction

    def inferior_conjunction_time(self, time, which='nearest'):
        """
        Time at inferior conjunction.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which mid-eclipse relative to the present ``time`` would you
            like to calculate. Default is nearest.
        Returns
        -------
        `~astropy.time.Time`
            Time of inferior conjunction
        """
        if not isinstance(time, Time):
            time = Time(time)

        phase_at_time = (abs(time - self.time_inferior_conjunction).to(u.day) %
                         self.period)
        next_eclipse = time + (self.period - phase_at_time)
        previous_eclipse = time - phase_at_time

        if which == 'next':
            return next_eclipse
        elif which == 'previous':
            return previous_eclipse
        elif which == 'nearest':
            if abs(next_eclipse - time) < abs(previous_eclipse - time):
                return next_eclipse
            else:
                return previous_eclipse
        else:
            raise ValueError('"which" kwarg must be "next", "previous" or '
                             '"nearest".')

    def ingress_time(self, time, which='nearest'):
        """
        Ingress time.

        Assumes that the difference between time of ingress and the time of
        mid-eclipse is equivalent to the difference between the time of egress
        to the time of egress (i.e. strictly true for circular orbits only).

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which ingress relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of ingress
        """
        if self.eclipse_duration is None:
            raise ValueError("An eclipse duration is required to compute "
                             "ingress times.")

        if not isinstance(time, Time):
            time = Time(time)

        phase_at_time = (abs(time - self.time_inferior_conjunction).to(u.day) %
                         self.period)
        next_ingress = time + (self.period - phase_at_time) - 0.5*self.eclipse_duration
        previous_ingress = time - phase_at_time - 0.5*self.eclipse_duration

        if which == 'next':
            return next_ingress
        elif which == 'previous':
            return previous_ingress
        elif which == 'nearest':
            if abs(next_ingress - time) < abs(previous_ingress - time):
                return next_ingress
            else:
                return previous_ingress
        else:
            raise ValueError('"which" kwarg must be "next", "previous" or '
                             '"nearest".')

    def egress_time(self, time, which='nearest'):
        """
        Egress time.

        Assumes that the difference between time of ingress and the time of
        mid-eclipse is equivalent to the difference between the time of egress
        to the time of egress (i.e. strictly true for circular orbits only).

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which egress relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of egress
        """
        if self.eclipse_duration is None:
            raise ValueError("An eclipse duration is required to compute "
                             "egress times.")

        if not isinstance(time, Time):
            time = Time(time)

        phase_at_time = (abs(time - self.time_inferior_conjunction).to(u.day) %
                         self.period)
        next_egress = time + (self.period - phase_at_time) + 0.5*self.eclipse_duration
        previous_egress = time - phase_at_time + 0.5*self.eclipse_duration

        if which == 'next':
            return next_egress
        elif which == 'previous':
            return previous_egress
        elif which == 'nearest':
            if abs(next_egress - time) < abs(previous_egress - time):
                return next_egress
            else:
                return previous_egress
        else:
            raise ValueError('"which" kwarg must be "next", "previous" or '
                             '"nearest".')

    def phase(self, time):
        """
        Orbital phase of secondary on [0, 1].

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time(s) of observations. This will be passed in as the first
            argument to the `~astropy.time.Time` initializer, so it can be
            anything that `~astropy.time.Time` will accept (including a
            `~astropy.time.Time` object).

        Returns
        -------
        orbital_phase : float or array
            Orbital phase of secondary at ``times``.
        """
        if not isinstance(time, Time):
            time = Time(time)

        orbital_phase = (abs(time - self.time_inferior_conjunction).to(u.day) /
                         self.period) % 1
        return orbital_phase

class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """
