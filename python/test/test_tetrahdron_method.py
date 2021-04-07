import numpy as np
from kspclib import (get_thm_relative_grid_addresses,
                     get_all_grid_addresses,
                     get_all_grgrid_addresses,
                     get_double_grid_index,
                     get_double_grid_address,
                     get_double_grgrid_index,
                     get_double_grgrid_address,
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


def test_get_thm_integration_weight_nacl(nacl_lattice,
                                         nacl_phonon_frequences_101010):
    dos_str_df_1 = """0 6.695122056070627165e-05 8.259314625201316929e-07
1 8.911857347254222017e-02 2.522323511872215027e-02
2 4.685384246034071665e-01 2.630086572378442233e-01
3 1.697958549809626128e+00 1.353336182540944232e+00
4 2.388574749472669012e+00 2.587628545445126438e+00
5 3.892363708541792366e+00 4.944242466208688569e+00
6 5.985385823393145621e-01 5.700405664774345738e+00
7 6.634899679165588704e-02 5.991230128112197129e+00"""
    dos_ref = np.fromstring(dos_str_df_1, sep=' ').reshape(-1, 3)
    freqs = nacl_phonon_frequences_101010
    mesh = [10, 10, 10]
    num_gps = np.prod(mesh)
    num_bands = 6
    shift = [0, 0, 0]
    rec_lat = np.linalg.inv(nacl_lattice)
    grid_addresses = get_all_grid_addresses(mesh)
    relative_addresses = get_thm_relative_grid_addresses(rec_lat)
    fpoints = dos_ref[:, 0]
    dos = np.zeros_like(fpoints)
    acc = np.zeros_like(fpoints)
    for ga in grid_addresses:
        tetrahedra_gps = _get_tetrahedra_grid_indices(
            ga + relative_addresses, mesh, shift)
        tetrahedra_freqs = freqs[tetrahedra_gps]
        for i, fpt in enumerate(fpoints):
            for j in range(num_bands):
                dos[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j])
                acc[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j], function='J')

    dos_dat = np.array([fpoints, dos / num_gps, acc / num_gps]).T
    # np.savetxt('dos.dat', dos_dat)
    np.testing.assert_allclose(dos_dat, dos_ref, atol=1e-5)


def test_get_thm_relative_grgrid_addresses_nacl():
    P = [[0, 0, 1],
         [0, -1, -1],
         [1, 1, 0]]
    grid_matrix = [[2, -2, 2],
                   [2, 2, -2],
                   [-2, 2, 2]]
    # Non-spglib-conventional choice of primitive basis vectors
    # in row vectors.
    nacl_lattice = np.array([[0.5, 0.5, 0],
                             [0, 0.5, 0.5],
                             [0.5, 0, 0.5]]) * 4.0729350500
    rec_lat = np.linalg.inv(nacl_lattice)
    microzone = np.dot(rec_lat, np.linalg.inv(grid_matrix))
    relative_addresses = get_thm_relative_grid_addresses(microzone)
    gr_relative_addresses = np.dot(relative_addresses, np.transpose(P))

    t1 = [[0, 0, 0], [0, 0, 1], [0, -1, 2], [1, -2, 2]]
    t24 = [[0, 0, 0], [-1, 2, -2], [0, 1, -2], [0, 0, -1]]
    np.testing.assert_array_equal(gr_relative_addresses[0], t1)
    np.testing.assert_array_equal(gr_relative_addresses[23], t24)


def test_get_thm_integration_weight_tio2(tio2_lattice,
                                         tio2_phonon_frequences_grg):

    dos_str_df_1 = """2.5 0.33432705 0.22741990
5.0 2.51042492 2.93972973
7.5 0.85214354 5.31944265
10.0 1.22209820 8.51003829
12.5 1.36331347 10.72008644
15.0 1.84840716 13.72386364
17.5 0.25849856 14.97997129
20.0 0.05704757 16.04711651
22.5 0.40515558 16.44268108
25.0 0.00814283 17.99943050"""
    dos_ref = np.fromstring(dos_str_df_1, sep=' ').reshape(-1, 3)
    freqs = tio2_phonon_frequences_grg
    grid_matrix = [[0, 11, 11],
                   [11, 0, 11],
                   [4, 4, 0]]
    D_diag = [1, 11, 88]
    P = [[0, -1, 3],
         [1, 0, 0],
         [-4, 4, -11]]
    num_gps = np.prod(D_diag)
    num_bands = freqs.shape[1]
    shift = [0, 0, 0]
    rec_lat = np.linalg.inv(tio2_lattice)
    microzone = np.dot(rec_lat, np.linalg.inv(grid_matrix))
    grid_addresses = get_all_grgrid_addresses(D_diag)
    # print(grid_addresses[:20])
    relative_addresses = get_thm_relative_grid_addresses(microzone)
    gr_relative_addresses = np.dot(relative_addresses, np.transpose(P))
    fpoints = dos_ref[:, 0]
    dos = np.zeros_like(fpoints)
    acc = np.zeros_like(fpoints)
    for ga in grid_addresses:
        tetrahedra_gps = _get_tetrahedra_grgrid_indices(
            ga + gr_relative_addresses, D_diag, shift)
        tetrahedra_freqs = freqs[tetrahedra_gps]
        for i, fpt in enumerate(fpoints):
            for j in range(num_bands):
                dos[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j])
                acc[i] += get_thm_integration_weight(
                    fpt, tetrahedra_freqs[:, :, j], function='J')

    dos_dat = np.array([fpoints, dos / num_gps, acc / num_gps]).T
    # for line in dos_dat:
    #     print("%.1f %.8f %.8f" % tuple(line))
    np.testing.assert_allclose(dos_dat, dos_ref, atol=1e-5)


def _get_frequency_points(freqs, df=0.1):
    fmax = np.amax(freqs)
    fmin = 0.0429161884
    fpoints = np.arange(fmin, fmax, df)
    return fpoints


def _get_tetrahedra_grid_indices(tetrahedra_ga, mesh, shift):
    tetrahedra_gps = np.zeros((24, 4), dtype='uintp', order='C')
    for j in range(24):
        for k in range(4):
            ga_d = get_double_grid_address(
                tetrahedra_ga[j, k], mesh, shift=shift)
            tetrahedra_gps[j, k] = get_double_grid_index(ga_d, mesh)
    return tetrahedra_gps


def _get_tetrahedra_grgrid_indices(tetrahedra_ga, D_diag, PS):
    tetrahedra_gps = np.zeros((24, 4), dtype='int_', order='C')
    for j in range(24):
        for k in range(4):
            ga_d = get_double_grgrid_address(
                tetrahedra_ga[j, k], D_diag, PS=PS)
            tetrahedra_gps[j, k] = get_double_grgrid_index(ga_d, D_diag, PS=PS)
    return tetrahedra_gps
