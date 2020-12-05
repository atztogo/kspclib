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
    """row vectors"""
    lattice = [[ -1.888070425000000, 1.888070425000000, 4.790243149999999 ],
               [ 1.888070425000000, -1.888070424999999, 4.790243149999999 ],
               [ 1.888070425000000, 1.888070424999999, -4.790243149999999 ]]
    return lattice
