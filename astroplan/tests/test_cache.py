from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..cache import pickle_data, unpickle_data

def test_pickling():
    pickle_this = dict(a=1, b=2, c=3)
    file_name = "test.pkl"
    pickle_data(pickle_this, file_name)
    unpickled_data = unpickle_data(file_name)
    assert unpickled_data == pickle_this

