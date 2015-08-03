# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Manage cached files.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.config.paths import get_cache_dir
import os
import sys
from astropy.extern.six.moves import cPickle as pickle

__all__ = ["unpickle_data", "pickle_data"]

def _find_or_create_astroplan_dir(dirnm):
    """
    Adapted for astroplan from astropy.config.paths._find_or_create_astropy_dir
    """
    innerdir = os.path.join(get_cache_dir(), 'astroplan')
    maindir = os.path.join(get_cache_dir(), 'astroplan', dirnm)

    if not os.path.exists(maindir):
        # first create .astropy/astroplan dir if needed
        if not os.path.exists(innerdir):
            try:
                os.mkdir(innerdir)
            except OSError:
                if not os.path.isdir(innerdir):
                    raise
        elif not os.path.isdir(innerdir):
            msg = 'Intended astroplan directory {0} is actually a file.'
            raise IOError(msg.format(innerdir))

        try:
            os.mkdir(maindir)
        except OSError:
            if not os.path.isdir(maindir):
                raise

    elif not os.path.isdir(maindir):
        msg = 'Intended astroplan {0} directory {1} is actually a file.'
        raise IOError(msg.format(dirnm, maindir))

    return os.path.abspath(maindir)

def pickle_data(data, file_name):
    """
    Pickle the object ``data`` into the astroplan cache, with the name
    ``file_name``.

    Parameters
    ----------
    data : obj
        Thing to pickle

    file_name : str
        Name of the pickle file (no dir needed)
    """
    pickle_path = os.path.join(_find_or_create_astroplan_dir('data'), file_name)
    with open(pickle_path, 'wb') as f:
        pickle.dump(data, f)

def unpickle_data(file_name):
    """
    Unpickle the object stored in the astroplan cache with the name
    ``file_name``.

    Parameters
    ----------
    file_name : str
        Name of the pickle file (no dir needed)

    Returns
    -------
    data : obj
        Stored file.
    """
    pickle_path = os.path.join(_find_or_create_astroplan_dir('data'), file_name)
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
    return data
