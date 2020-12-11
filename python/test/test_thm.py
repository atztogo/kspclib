import numpy as np
from kspclib import (get_thm_relative_grid_addresses,
                     get_all_grid_addresses,
                     get_double_grid_point,
                     get_double_grid_address,
                     get_thm_integration_weight)


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


def test_get_thm_integration_weight(nacl_lattice,
                                    nacl_phonon_frequences_101010):
    dos_str_df_1 = """0 6.695122056070627165e-05 8.259314625201316929e-07
1 8.911857347254222017e-02 2.553382929530613465e-02
2 4.685384246034071665e-01 2.649562989742248464e-01
3 1.697958549809626128e+00 1.362664683321239689e+00
4 2.388574749472669012e+00 2.596688544653257491e+00
5 3.892363708541792366e+00 4.946293094871299090e+00
6 5.985385823393145621e-01 5.703047775018656118e+00
7 6.634899679165588704e-02 5.991502313827999693e+00"""

    freqs = nacl_phonon_frequences_101010
    mesh = [10, 10, 10]
    num_gps = np.prod(mesh)
    num_bands = 6
    shift = [0, 0, 0]
    df = 1.0
    rec_lat = np.linalg.inv(nacl_lattice)
    grid_addresses = get_all_grid_addresses(mesh)
    relative_addresses = get_thm_relative_grid_addresses(rec_lat)
    fpoints = _get_frequency_points(freqs, df=df)
    dos = np.zeros_like(fpoints)
    acc = np.zeros_like(fpoints)
    for ga in grid_addresses:
        tetrahedra_gps = _get_tetrahedra_grid_points(
            ga + relative_addresses, mesh, shift)
        tetrahedra_freqs = freqs[tetrahedra_gps]
        for i, fpt in enumerate(fpoints):
            for j in range(num_bands):
                dos[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j])
                acc[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j], function='J')

    dos_dat = np.array([fpoints, dos / num_gps, acc / num_gps]).T
    dos_ref = np.fromstring(dos_str_df_1, sep=' ')
    np.testing.assert_allclose(dos_dat.ravel(), dos_ref, atol=1e-5)
    # np.savetxt('dos.dat', dos_dat)


def _get_frequency_points(freqs, df=0.1):
    fmax = np.amax(freqs)
    fmin = 0
    fpoints = np.arange(fmin, fmax, df)
    return fpoints


def _get_tetrahedra_grid_points(tetrahedra_ga, mesh, shift):
    tetrahedra_gps = np.zeros((24, 4), dtype='uintp', order='C')
    for j in range(24):
        for k in range(4):
            ga_d = get_double_grid_address(
                tetrahedra_ga[j, k], mesh, shift=shift)
            tetrahedra_gps[j, k] = get_double_grid_point(ga_d, mesh)
    return tetrahedra_gps
