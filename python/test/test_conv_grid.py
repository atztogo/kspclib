import numpy as np
from kspclib import (get_all_grid_addresses, get_double_grid_point,
                     get_double_grid_address)

addresses_234 = np.array([[0, 0, 0],
                          [1, 0, 0],
                          [0, 1, 0],
                          [1, 1, 0],
                          [0, -1, 0],
                          [1, -1, 0],
                          [0, 0, 1],
                          [1, 0, 1],
                          [0, 1, 1],
                          [1, 1, 1],
                          [0, -1, 1],
                          [1, -1, 1],
                          [0, 0, 2],
                          [1, 0, 2],
                          [0, 1, 2],
                          [1, 1, 2],
                          [0, -1, 2],
                          [1, -1, 2],
                          [0, 0, -1],
                          [1, 0, -1],
                          [0, 1, -1],
                          [1, 1, -1],
                          [0, -1, -1],
                          [1, -1, -1]])


def test_get_all_grid_addresses():
    mesh = [2, 3, 4]
    grid_addresses = get_all_grid_addresses(mesh)
    np.testing.assert_array_equal(grid_addresses, addresses_234)


def test_get_double_grid_point():
    mesh = [2, 3, 4]
    for shift in np.ndindex((2, 2, 2)):
        for i, address in enumerate(addresses_234):
            gp = get_double_grid_point(address * 2 + shift, mesh)
            assert i == gp


def test_get_double_grid_address():
    mesh = [2, 3, 4]
    for shift in np.ndindex((2, 2, 2)):
        for i, address in enumerate(addresses_234):
            gp = get_double_grid_point(address * 2 + shift, mesh)
            ga = get_double_grid_address(address, mesh, shift=shift)
            ga_gp = get_double_grid_point(ga, mesh)
            assert gp == ga_gp
