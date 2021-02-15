import numpy as np
from kspclib import get_reciprocal_point_group


def test_tipn3_rotations(tipn3_direct_rotations, tipn3_reciprocal_rotations):
    rec_rots = get_reciprocal_point_group(tipn3_direct_rotations)
    np.testing.assert_array_equal(rec_rots, tipn3_reciprocal_rotations)


def test_tipn3_rotations_no_time_reversal(tipn3_direct_rotations,
                                          tipn3_reciprocal_rotations):
    rec_rots = get_reciprocal_point_group(tipn3_direct_rotations,
                                          is_time_reversal=False)
    np.testing.assert_array_equal(np.vstack([rec_rots, -1 * rec_rots]),
                                  tipn3_reciprocal_rotations)


def test_tio2_rotations(tio2_direct_rotations, tio2_reciprocal_rotations):
    rec_rots = get_reciprocal_point_group(tio2_direct_rotations)
    np.testing.assert_array_equal(rec_rots, tio2_reciprocal_rotations)


def test_tio2_rotations_nonprim(tio2_direct_rotations,
                                tio2_reciprocal_rotations):
    rots = tio2_direct_rotations
    rots_nonprim = np.vstack([rots, rots])
    rec_rots = get_reciprocal_point_group(rots_nonprim)
    np.testing.assert_array_equal(rec_rots, tio2_reciprocal_rotations)
