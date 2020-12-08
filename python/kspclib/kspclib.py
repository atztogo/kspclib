# Copyright (C) 2020 Atsushi Togo
# All rights reserved.
#
# This file is part of kspclib.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the kspclib project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

from . import _kspclib as ksp
import numpy as np


def get_version():
    return tuple(ksp.version())


def get_all_grid_addresses(mesh):
    """Return all grid addresses for mesh

    Parameters
    ----------
    mesh : array_like
        Conventional regular mesh for grid sampling.
        shape=(3,), dtype='intc'

    Returns
    -------
    grid_address : ndarray
        Grid addresses of all grid points corresponding to input mesh.
        shape=(all_grid_points, 3), dtype='intc'

    """

    grid_address = np.zeros((np.prod(mesh), 3), dtype='intc', order='C')
    ksp.all_grid_addresses(grid_address, np.array(mesh, dtype='intc'))
    return grid_address


def get_grid_point_double_mesh(address_double, mesh):
    """Return grid point index of grid address of mesh

    Parameters
    ----------
    address_double : array_like
        Grid address in double grid method.
        shape=(3,), dtype='intc'
    mesh : array_like
        Conventional regular mesh for grid sampling.
        shape=(3,), dtype='intc'

    Returns
    -------
    grid_point : int
        Grid point index.

    """
    return ksp.grid_point_double_mesh(np.array(address_double, dtype='intc'),
                                      np.array(mesh, dtype='intc'))


def get_grid_address_double_mesh(address, mesh, is_shift=None):
    """Return grid point index of grid address of mesh

    Parameters
    ----------
    address : array_like
        Grid address.
        shape=(3,), dtype='intc'
    mesh : array_like
        Conventional regular mesh for grid sampling.
        shape=(3,), dtype='intc'
    is_shift : array_like, optional
        Half grid shift for conventional regular mesh along reciprocal basis
        vector directions. 0 and 1 mean no shift and half shift, recpectively.
        shape=(3,), dtype='intc'

    Returns
    -------
    address_double : ndarray
        Grid address in double grid method.
        shape=(3,), dtype='intc'

    """

    address_double = np.zeros(3, dtype='intc', order='C')
    if is_shift is None:
        _is_shift = np.zeros(3, dtype='intc', order='C')
    else:
        _is_shift = np.array(is_shift, dtype='intc')
    ksp.grid_address_double_mesh(address_double,
                                 np.array(address, dtype='intc'),
                                 np.array(mesh, dtype='intc'),
                                 _is_shift)
    return address_double


def get_thm_relative_grid_addresses(rec_lattice):
    """Return relative grid addresses of 24 tetrahedra

    Parameters
    ----------
    rec_lattice : array_like
       Reciprocal basis vectors in column vectors.
       shape=(3, 3), dtype='double', order='C'

    Returns
    -------
    relative_addresses : ndarray
       Grid address shifts corresponding to 24 tetrahedra surrounding
       a grid point for conventional regular grid.
       shape=(24, 4, 3), dtype='intc', order='C'

    """

    relative_addresses = np.zeros((24, 4, 3), dtype='intc', order='C')
    ksp.thm_relative_grid_addresses(
        relative_addresses,
        np.array(rec_lattice, dtype='double', order='C'))
    return relative_addresses


def get_thm_integration_weight(omega, tetrahedra_omegas, function='I'):
    """Return tetheradron method integration weight for a grid point

    Parameters
    ----------
    omega : float
        Energy where integration weight is computed.
    tetrahedra_omegas : array_like
        Energies of four vertices of 24 tetrahedra. These energies are those
        at the grid points as given by ``get_thm_relative_grid_addresses``.
        shape=(24, 4), dtype='double', order='C'
    function : str, optional
        'I' for delta function and 'J' for Heaviside function. Default is 'I'.

    Returns
    -------
    integration_weight : float
        Integration weight for a grid point.

    """

    iw = ksp.thm_integration_weight(
        float(omega),
        np.array(tetrahedra_omegas, dtype='double', order='C'),
        str(function.upper()))
    return iw


def get_snf3x3(A):
    """Return Smith normal form of 3x3 integer matrix

    Parameters
    ----------
    A : array_like
        Integer transformation matrix from basis vectors of microzone to
        those of primitive basis vectors.
        shape=(3, 3), dtype='int_', order='C'

    returns
    -------
    snf : dict
        D, P, Q of Smith normal form of 3x3 integer matrix. The dict keys are
        'D', 'P', 'Q', respectively.
        D, P, Q: shape=(3, 3), dtype='int_', order='C'

    """

    DPQ = np.zeros((3, 3, 3), dtype='int_', order='C')
    succeeded = ksp.snf3x3(DPQ, np.array(A, dtype='int_', order='C'))

    if succeeded:
        return {'D': np.array(DPQ[0], dtype='int_', order='C'),
                'P': np.array(DPQ[1], dtype='int_', order='C'),
                'Q': np.array(DPQ[2], dtype='int_', order='C')}
    else:
        return None


def sanity_check_rotations(rotations, grid_matrix=None, D=None, Q=None):
    """Check compatibility of grid generation matrix against rotations

    Parameters
    ----------
    D, Q : array_like, optional
        D and Q of Smith normal form of grid matrix. Default is None.
        shape=(3, 3), dtype='int_', order='C'.
        For D, a sequence of diagonal elements with shape=(3,) is accepted.
    grid_matrix : array_like, optional
        Grid generation matrix. Default is None.
        shape=(3, 3), dtype='int_', order='C'
    rotations : array_like
        Reciprocal rotation matrices.
        shape=(num_rot, 3, 3), dtype='intc', order='C'

    returns
    -------
    bool
        True: compatible, False: incompatible.

    """

    if grid_matrix is not None:
        snf = get_snf3x3(grid_matrix)
        _D = snf['D']
        _Q = snf['Q']
    elif D is not None and Q is not None:
        if len(np.ravel(D)) == 3:
            _D = np.diag(np.ravel(D))
        else:
            _D = D
        _Q = Q
    else:
        msg = "grid_matrix or D and Q unspecified."
        raise RuntimeError(msg)

    is_compatible = ksp.sanity_check_rotations(
        np.array(_D, dtype='int_', order='C'),
        np.array(_Q, dtype='int_', order='C'),
        np.array(rotations, dtype='intc', order='C'))

    return is_compatible
