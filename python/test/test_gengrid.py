import numpy as np
from kspclib import (get_snf3x3, )


def test_get_snf3x3():
    A = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype='int_') * 100
    snf = get_snf3x3(A)
    _D = np.dot(snf['P'], np.dot(A, snf['Q']))
    np.testing.assert_array_equal(snf['D'], _D)
