import numpy as np
from kspclib import (get_thm_relative_grid_addresses, )


def test_get_thm_relative_grid_addresses_nacl(nacl_lattice):
    rec_lat = np.linalg.inv(nacl_lattice)
    relative_addresses = get_thm_relative_grid_addresses(rec_lat)
    t1 = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]
    t24 = [[0, 0, 0], [-1, -1, -1], [-1, -1, 0], [-1, 0, 0]]
    np.testing.assert_array_equal(relative_addresses[0], t1)
    np.testing.assert_array_equal(relative_addresses[23], t24)


def test_get_thm_relative_grid_addresses_sio2(sio2_lattice):
    rec_lat = np.linalg.inv(sio2_lattice)
    relative_addresses = get_thm_relative_grid_addresses(rec_lat)
    t1 = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]
    t24 = [[0, 0, 0], [-1, -1, -1], [-1, -1, 0], [-1, 0, 0]]
    np.testing.assert_array_equal(relative_addresses[0], t1)
    np.testing.assert_array_equal(relative_addresses[23], t24)


def test_get_thm_relative_grid_addresses_tio2(tio2_lattice):
    rec_lat = np.linalg.inv(tio2_lattice)
    relative_addresses = get_thm_relative_grid_addresses(rec_lat)
    t1 = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1]]
    t24 = [[0, 0, 0], [0, 0, -1], [-1, -1, 0], [-1, 0, 0]]
    np.testing.assert_array_equal(relative_addresses[0], t1)
    np.testing.assert_array_equal(relative_addresses[23], t24)
