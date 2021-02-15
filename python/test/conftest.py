import os
import pytest
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope='session')
def nacl_lattice():
    x = 5.6903014761756712 / 2
    lattice = [[0, x, x], [x, 0, x], [x, x, 0]]
    return lattice


@pytest.fixture(scope='session')
def sio2_lattice():
    lattice = [[4.65, 0, 0],
               [0, 4.75, 0],
               [0, 0, 3.25]]
    return lattice


@pytest.fixture(scope='session')
def tio2_lattice():
    """Primitive cell basis vectors of TiO2 anatase in row vectors"""
    lattice = [[-1.888070425000000, 1.888070425000000, 4.790243149999999],
               [1.888070425000000, -1.888070424999999, 4.790243149999999],
               [1.888070425000000, 1.888070424999999, -4.790243149999999]]
    return lattice


@pytest.fixture(scope='session')
def nacl_phonon_frequences_101010():
    """NaCl phonon frequency data

    Phonons are sampled on 10x10x10 regular grid without shift.
    Number of bands is six.

    """
    filename = os.path.join(current_dir, "NaCl-freqs-101010.dat")
    freqs = np.loadtxt(filename).reshape(-1, 6)
    return freqs


@pytest.fixture(scope='session')
def tio2_phonon_frequences_grg():
    """TiO2 phonon frequency data

    Phonons are sampled on the grid matrix defined below.

    grid_matrix = [[0, 11, 11],
                   [11, 0, 11],
                   [4, 4, 0]]

    """
    filename = os.path.join(current_dir, "TiO2-freqs-grg.dat")
    freqs = np.loadtxt(filename).reshape(-1, 18)
    return freqs


@pytest.fixture(scope='session')
def tio2_direct_rotations():
    rotations = np.reshape([1, 0, 0, 0, 1, 0, 0, 0, 1,
                            0, 1, 0, 0, 1, -1, -1, 1, 0,
                            0, 1, -1, 1, 0, -1, 0, 0, -1,
                            1, 0, -1, 1, 0, 0, 1, -1, 0,
                            -1, 0, 0, -1, 0, 1, -1, 1, 0,
                            0, -1, 0, -1, 0, 0, 0, 0, -1,
                            0, -1, 1, 0, -1, 0, 1, -1, 0,
                            -1, 0, 1, 0, -1, 1, 0, 0, 1,
                            -1, 0, 0, 0, -1, 0, 0, 0, -1,
                            0, -1, 0, 0, -1, 1, 1, -1, 0,
                            0, -1, 1, -1, 0, 1, 0, 0, 1,
                            -1, 0, 1, -1, 0, 0, -1, 1, 0,
                            1, 0, 0, 1, 0, -1, 1, -1, 0,
                            0, 1, 0, 1, 0, 0, 0, 0, 1,
                            0, 1, -1, 0, 1, 0, -1, 1, 0,
                            1, 0, -1, 0, 1, -1, 0, 0, -1], (-1, 3, 3))
    return rotations


@pytest.fixture(scope='session')
def tio2_reciprocal_rotations():
    rotations = np.reshape([1, 0, 0, 0, 1, 0, 0, 0, 1,
                            0, 0, -1, 1, 1, 1, 0, -1, 0,
                            0, 1, 0, 1, 0, 0, -1, -1, -1,
                            1, 1, 1, 0, 0, -1, -1, 0, 0,
                            -1, -1, -1, 0, 0, 1, 0, 1, 0,
                            0, -1, 0, -1, 0, 0, 0, 0, -1,
                            0, 0, 1, -1, -1, -1, 1, 0, 0,
                            -1, 0, 0, 0, -1, 0, 1, 1, 1,
                            -1, 0, 0, 0, -1, 0, 0, 0, -1,
                            0, 0, 1, -1, -1, -1, 0, 1, 0,
                            0, -1, 0, -1, 0, 0, 1, 1, 1,
                            -1, -1, -1, 0, 0, 1, 1, 0, 0,
                            1, 1, 1, 0, 0, -1, 0, -1, 0,
                            0, 1, 0, 1, 0, 0, 0, 0, 1,
                            0, 0, -1, 1, 1, 1, -1, 0, 0,
                            1, 0, 0, 0, 1, 0, -1, -1, -1], (-1, 3, 3))
    return rotations


@pytest.fixture(scope='session')
def tipn3_direct_rotations():
    rotations = np.reshape([1, 0, 0, 0, 1, 0, 0, 0, 1,
                            -1, 0, 0, 0, 0, 1, 0, 1, 0,
                            -1, 0, 0, 0, 1, 0, 0, 0, 1,
                            1, 0, 0, 0, 0, 1, 0, 1, 0], (-1, 3, 3))
    return rotations


@pytest.fixture(scope='session')
def tipn3_reciprocal_rotations():
    rotations = np.reshape([1, 0, 0, 0, 1, 0, 0, 0, 1,
                            -1, 0, 0, 0, 0, 1, 0, 1, 0,
                            -1, 0, 0, 0, 1, 0, 0, 0, 1,
                            1, 0, 0, 0, 0, 1, 0, 1, 0,
                            -1, 0, 0, 0, -1, 0, 0, 0, -1,
                            1, 0, 0, 0, 0, -1, 0, -1, 0,
                            1, 0, 0, 0, -1, 0, 0, 0, -1,
                            -1, 0, 0, 0, 0, -1, 0, -1, 0], (-1, 3, 3))
    return rotations
