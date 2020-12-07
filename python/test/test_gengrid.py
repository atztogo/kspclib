import numpy as np
from kspclib import (get_snf3x3, )


def test_get_snf3x3():
    A = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype='int_') * 100
    D, P, Q = get_snf3x3(A)
    _D = np.dot(P, np.dot(A, Q))
    np.testing.assert_array_equal(D, _D)
