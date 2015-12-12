import cython
import logging
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "spglib_interface.cpp":
        void get_reciprocal_mesh_interface(int *, double *, double *, int *, unsigned int, int *, int *, int *, int, double)
        
@cython.boundscheck(False)
@cython.wraparound(False)

def get_reciprocal_mesh(np.ndarray[int, ndim=1, mode="c"] mesh not None, np.ndarray[double, ndim=2, mode="c"] lattice not None, np.ndarray[double, ndim=2, mode="c"] positions not None, np.ndarray[int, ndim=1, mode="c"] anumbers not None, np.ndarray[int, ndim=1, mode="c"] is_shifted not None, np.ndarray[int, ndim=2, mode="c"] mesh_points not None, np.ndarray[int,ndim=1,mode="c"] mapping not None, is_time_reversal=True, symprec=1e-5):
        num_atoms=anumbers.shape[0]
        get_reciprocal_mesh_interface(&mesh[0], &lattice[0,0], &positions[0,0], &anumbers[0], num_atoms, &is_shifted[0], &mesh_points[0,0], &mapping[0], is_time_reversal, symprec)
