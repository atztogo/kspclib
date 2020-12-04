from kspclib import get_all_grid_addresses


def test_get_all_grid_addresses():
    mesh = [2, 2, 2]
    grid_addresses = get_all_grid_addresses(mesh)
    print(grid_addresses)
