import numpy as np
from kspclib import (get_snf3x3, snf_transform_rotations,
                     get_all_grgrid_addresses)

# (16, 3, 3)
tio2_rots = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 1],
    [0, -1, 0],
    [0, 1, 0],
    [1, 0, 0],
    [-1, -1, -1],
    [1, 1, 1],
    [0, 0, -1],
    [-1, 0, 0],
    [-1, -1, -1],
    [0, 0, 1],
    [0, 1, 0],
    [0, -1, 0],
    [-1, 0, 0],
    [0, 0, -1],
    [0, 0, 1],
    [-1, -1, -1],
    [1, 0, 0],
    [-1, 0, 0],
    [0, -1, 0],
    [1, 1, 1],
    [-1, 0, 0],
    [0, -1, 0],
    [0, 0, -1],
    [0, 0, 1],
    [-1, -1, -1],
    [0, 1, 0],
    [0, -1, 0],
    [-1, 0, 0],
    [1, 1, 1],
    [-1, -1, -1],
    [0, 0, 1],
    [1, 0, 0],
    [1, 1, 1],
    [0, 0, -1],
    [0, -1, 0],
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 1],
    [-1, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [-1, -1, -1]]
tio2_grid_mat = [[0, 5, 5],
                 [5, 0, 5],
                 [2, 2, 0]]
tio2_snf = {'D_diag': [1, 5, 20],
            'P': [[0, 1, -2],
                  [1, 0, 0],
                  [-2, 2, -5]],
            'Q': [[1, -5, -9],
                  [0, 0, -1],
                  [0, 1, 1]]}

# (16, 3, 3)
tio2_transformed_rots = [  # by phonopy's GeneralizedRegularGridPoints
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [-4, 3, 2],
    [5, -4, -2],
    [-20, 16, 9],
    [-9, 8, 4],
    [0, -1, 0],
    [-20, 20, 9],
    [-4, 5, 2],
    [-5, 4, 2],
    [0, 4, 1],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -4, -1],
    [4, -5, -2],
    [-5, 4, 2],
    [20, -20, -9],
    [9, -8, -4],
    [0, -1, 0],
    [20, -16, -9],
    [4, -3, -2],
    [5, -4, -2],
    [0, 0, -1],
    [-1, 0, 0],
    [0, -1, 0],
    [0, 0, -1],
    [4, -3, -2],
    [-5, 4, 2],
    [20, -16, -9],
    [9, -8, -4],
    [0, 1, 0],
    [20, -20, -9],
    [4, -5, -2],
    [5, -4, -2],
    [0, -4, -1],
    [1, 0, 0],
    [0, -1, 0],
    [0, 4, 1],
    [-4, 5, 2],
    [5, -4, -2],
    [-20, 20, 9],
    [-9, 8, 4],
    [0, 1, 0],
    [-20, 16, 9],
    [-4, 3, 2],
    [-5, 4, 2],
    [0, 0, 1]]

addresses_244 = np.array([[0, 0, 0],
                          [1, 0, 0],
                          [0, 1, 0],
                          [1, 1, 0],
                          [0, 2, 0],
                          [1, 2, 0],
                          [0, 3, 0],
                          [1, 3, 0],
                          [0, 0, 1],
                          [1, 0, 1],
                          [0, 1, 1],
                          [1, 1, 1],
                          [0, 2, 1],
                          [1, 2, 1],
                          [0, 3, 1],
                          [1, 3, 1],
                          [0, 0, 2],
                          [1, 0, 2],
                          [0, 1, 2],
                          [1, 1, 2],
                          [0, 2, 2],
                          [1, 2, 2],
                          [0, 3, 2],
                          [1, 3, 2],
                          [0, 0, 3],
                          [1, 0, 3],
                          [0, 1, 3],
                          [1, 1, 3],
                          [0, 2, 3],
                          [1, 2, 3],
                          [0, 3, 3],
                          [1, 3, 3]])


def test_get_snf3x3():
    A = tio2_grid_mat
    snf = get_snf3x3(A)
    _D = np.dot(snf['P'], np.dot(A, snf['Q']))
    np.testing.assert_array_equal(snf['D'], _D)
    for symbol in ('D_diag', 'P', 'Q'):
        np.testing.assert_array_equal(tio2_snf[symbol], snf[symbol])
    np.testing.assert_array_equal(snf['D'], np.diag(snf['D_diag']))


def test_snf_transform_rotations():
    """TiO2 anatase I4_1/amd (141)

    Primitive cell of body centred tetragonal crystal requires generalized
    regular mesh to let grid preserve symmetry. The conventional regular mesh
    (Q=I) fails although |a*|=|b*|.

    """

    rotations = np.reshape(tio2_rots, (16, 3, 3))
    # tio2_grid_mat = np.eye(3, dtype=int) * 4
    try:
        transformed_rots_1 = snf_transform_rotations(
            rotations, grid_matrix=tio2_grid_mat)
    except RuntimeError:
        assert False

    np.testing.assert_array_equal(
        transformed_rots_1.reshape(-1, 3), tio2_transformed_rots)

    try:
        transformed_rots_2 = snf_transform_rotations(
            rotations, D_diag=tio2_snf['D_diag'], Q=tio2_snf['Q'])
    except RuntimeError:
        assert False

    np.testing.assert_array_equal(
        transformed_rots_2.reshape(-1, 3), tio2_transformed_rots)

    try:
        transformed_rots_3 = snf_transform_rotations(
            rotations, D=np.diag(tio2_snf['D_diag']), Q=tio2_snf['Q'])
    except RuntimeError:
        assert False

    np.testing.assert_array_equal(
        transformed_rots_3.reshape(-1, 3), tio2_transformed_rots)

    try:
        transformed_rots_4 = snf_transform_rotations(
            rotations, D_diag=[4, 4, 8], Q=np.eye(3, dtype=int))
    except RuntimeError:
        assert True


def test_get_all_grgrid_addresses():
    D_diag = [2, 4, 4]
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    # for v in grgrid_addresses:
    #     print("[%d, %d, %d]," % tuple(v))
    np.testing.assert_array_equal(grgrid_addresses, addresses_244)
