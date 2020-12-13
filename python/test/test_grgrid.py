import numpy as np
from kspclib import (get_snf3x3, snf_transform_rotations,
                     get_all_grgrid_addresses,
                     get_double_grgrid_address,
                     get_grgrid_index,
                     get_double_grgrid_index,
                     niggli_reduce)

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


# def test_get_all_grgrid_addresses_for_q(tio2_lattice):
#     D_diag = tio2_snf['D_diag']
#     Q = np.array(tio2_snf['Q'])
#     grgrid_addresses = get_all_grgrid_addresses(D_diag)
#     q = np.dot(grgrid_addresses / np.array(D_diag, dtype=float), Q.T)
#     q -= np.rint(q)
#     from phonopy.structure.atoms import PhonopyAtoms
#     from phonopy.interface.vasp import write_vasp
#     cell = PhonopyAtoms(cell=tio2_lattice,
#                         scaled_positions=q,
#                         numbers=[1, ] * len(q))
#     write_vasp("POSCAR", cell)


def test_rotate_all_grgrid_addresses():
    D_diag = tio2_snf['D_diag']
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    ids = np.arange(len(grgrid_addresses))
    for r in np.reshape(tio2_transformed_rots, (-1, 3, 3)):
        rot_addresses = np.dot(grgrid_addresses, r.T)
        gps = [get_grgrid_index(adrs, D_diag)
               for adrs in rot_addresses]
        np.testing.assert_array_equal(np.sort(gps), ids)


def test_rotate_all_grgrid_double_addresses():
    """

    Primitive cell of TiO2 anataze, body centred tetragonal.
    Shifts have to be limited to (0, 0, 0), (1, 1, 0), (0, 0, 1),
    (1, 1, 1). Other shifts can fail rotating grid points, but not
    necessarily always failing.

    """

    D_diag = tio2_snf['D_diag']
    P = tio2_snf['P']
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    ids = np.arange(len(grgrid_addresses))
    for shift in ((0, 0, 0), (1, 1, 0), (0, 0, 1), (1, 1, 1)):
        PS = np.dot(P, shift)
        d_addresses = [get_double_grgrid_address(adrs, D_diag, PS=PS)
                       for adrs in grgrid_addresses]
        for r in np.reshape(tio2_transformed_rots, (-1, 3, 3)):
            rot_addresses = np.dot(d_addresses, r.T)
            gps = [get_double_grgrid_index(adrs, D_diag, PS=PS)
                   for adrs in rot_addresses]
            np.testing.assert_array_equal(np.sort(gps), ids)


def test_get_double_grid_index():
    D_diag = tio2_snf['D_diag']
    P = tio2_snf['P']
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    for shift in np.ndindex((2, 2, 2)):
        PS = np.dot(P, shift)
        for i, address in enumerate(grgrid_addresses):
            d_ga = get_double_grgrid_address(address, D_diag, PS=PS)
            d_gp = get_double_grgrid_index(d_ga, D_diag, PS=PS)
            assert i == d_gp


def test_get_grid_index():
    D_diag = tio2_snf['D_diag']
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    for i, address in enumerate(grgrid_addresses):
        s_gp = get_grgrid_index(address, D_diag)
        assert i == s_gp


def test_get_double_grid_address():
    D_diag = tio2_snf['D_diag']
    P = tio2_snf['P']
    grgrid_addresses = get_all_grgrid_addresses(D_diag)
    for shift in np.ndindex((2, 2, 2)):
        PS = np.dot(P, shift)
        for i, address in enumerate(grgrid_addresses):
            gp = get_double_grgrid_index(address * 2 + PS, D_diag, PS=PS)
            ga = get_double_grgrid_address(address, D_diag, PS=PS)
            ga_gp = get_double_grgrid_index(ga, D_diag, PS=PS)
            assert gp == ga_gp
