import numpy as np
from kspclib import (get_all_grid_addresses, get_grid_point_double_mesh,
                     get_grid_address_double_mesh)

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


def test_get_grid_point_double_mesh():
    mesh = [2, 3, 4]
    for shift in np.ndindex((2, 2, 2)):
        for i, address in enumerate(addresses_234):
            gp = get_grid_point_double_mesh(address * 2 + shift, mesh)
            assert i == gp


def test_get_grid_address_double_mesh():
    mesh = [2, 3, 4]
    for shift in np.ndindex((2, 2, 2)):
        for i, address in enumerate(addresses_234):
            gp = get_grid_point_double_mesh(address * 2 + shift, mesh)
            ga = get_grid_address_double_mesh(address, mesh, shift)
            ga_gp = get_grid_point_double_mesh(ga, mesh)
            assert gp == ga_gp
