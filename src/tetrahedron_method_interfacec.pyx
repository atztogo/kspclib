import cython
import logging
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "tetrahedron_method_interface.cpp":
	void calc_dos_lin_tetra(double *, int *, int *, int *, double *, double *, int, int, int, double, double *, double *)

@cython.boundscheck(False)
@cython.wraparound(False)

def calc_density_of_states_interface(np.ndarray[double, ndim=2, mode="c"] energy not None, np.ndarray[int, ndim=2, mode="c"] grid_address not None, np.ndarray[int, ndim=1, mode="c"] grid_mapping_table not None, np.ndarray[int, ndim=1, mode="c"] mesh not None, np.ndarray[double, ndim=2, mode="c"] rec_basis not None, np.ndarray[double, ndim=1, mode="c"] energy_samples not None, int num_energy_samples, int num_bands, int num_kpoints_ibz, double volume, np.ndarray[double, ndim=2, mode="c"] dos not None, np.ndarray[double, ndim=2, mode="c"] int_dos not None):
	calc_dos_lin_tetra(&energy[0,0], &grid_address[0,0], &grid_mapping_table[0], &mesh[0], &rec_basis[0,0], &energy_samples[0], num_energy_samples, num_bands, num_kpoints_ibz, volume, &dos[0,0], &int_dos[0,0])
