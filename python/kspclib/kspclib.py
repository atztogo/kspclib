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
    """Return all single-grid addresses

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


def get_double_grid_address(address, mesh, shift=None):
    """Convert grid address plus shift to double-grid address

    address_double = 2 * address + shift, where values are reduced to
    be closet to 0, -mesh/2 < address_double <= mesh/2.

    Parameters
    ----------
    address : array_like
        Grid address.
        shape=(3,), dtype='intc'
    mesh : array_like
        Conventional regular mesh for grid sampling.
        shape=(3,), dtype='intc'
    shift : array_like, optional
        Half grid shift for conventional regular mesh along reciprocal basis
        vector directions. 0 and 1 mean no shift and half shift, recpectively.
        shape=(3,), dtype='intc'

    Returns
    -------
    address_double : ndarray
        Double-grid address.
        shape=(3,), dtype='intc'

    """

    address_double = np.zeros(3, dtype='intc', order='C')
    if shift is None:
        _shift = np.zeros(3, dtype='intc', order='C')
    else:
        _shift = np.array(shift, dtype='intc')
    ksp.double_grid_address(address_double,
                            np.array(address, dtype='intc'),
                            np.array(mesh, dtype='intc'),
                            _shift)
    return address_double


def get_double_grid_index(address_double, mesh):
    """Return grid point index of a double-grid address

    Parameters
    ----------
    address_double : array_like
        Double-grid address.
        shape=(3,), dtype='intc'
    mesh : array_like
        Conventional regular mesh for grid sampling.
        shape=(3,), dtype='intc'

    Returns
    -------
    grid_index : int
        Grid point index.

    """
    return ksp.double_grid_index(np.array(address_double, dtype='intc'),
                                 np.array(mesh, dtype='intc'))


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
        D, P, Q of Smith normal form of 3x3 integer matrix.
        The dict keys are ``D``, ``D_diag``, ``P``, ``Q``, respectively.
        D, P, Q : shape=(3, 3), dtype='int_', order='C'
        D_diag : Diagonal elements of D, shape=(3,), dtype='int_', order='C'

    """

    D_diag = np.zeros(3, dtype='int_', order='C')
    P = np.zeros((3, 3), dtype='int_', order='C')
    Q = np.zeros((3, 3), dtype='int_', order='C')
    succeeded = ksp.snf3x3(D_diag, P, Q, np.array(A, dtype='int_', order='C'))
    D = np.array(np.diag(D_diag), dtype='int_', order='C')

    if succeeded:
        return {'D_diag': D_diag, 'D': D, 'P': P, 'Q': Q}
    else:
        return None


def snf_transform_rotations(rotations,
                            grid_matrix=None,
                            D=None,
                            D_diag=None,
                            Q=None):
    """Transform rotations by SNF of grid generation matrix

    Reciprocal rotation matrices of usual reciprocal basis vectors (R) are
    transformed to those of reciprocal basis vectors transformed by D and Q
    of Smith normal form of grid generation matrix. The formula implemented is

        DQ^{-1}RQD^{-1}.

    Grid generation matrix has to be compatible with R. Unless satisfied,
    exception is raised.

    Parameters
    ----------
    D, D_diag, Q : array_like, optional
        D or diagonal elemets of D and Q of Smith normal form of grid matrix.
        Default is None.
        D, Q : shape=(3, 3), dtype='int_', order='C'.
        D_diag : shape=(3,), dtype='int_', order='C'.
    grid_matrix : array_like, optional
        Grid generation matrix. Default is None.
        shape=(3, 3), dtype='int_', order='C'
    rotations : array_like
        Reciprocal rotation matrices of usual reciprocal basis vectors.
        shape=(num_rot, 3, 3), dtype='intc', order='C'

    returns
    -------
    transformed_rotations : ndarray
        Transformed reciprocal rotation matrices.
        shape=(num_rot, 3, 3), dtype='int_', order='C'

    """

    if grid_matrix is not None:
        snf = get_snf3x3(grid_matrix)
        _D_diag = snf['D_diag']
        _Q = snf['Q']
    elif D_diag is not None and Q is not None:
        _D_diag = D_diag
        _Q = Q
    elif D is not None and Q is not None:
        _D_diag = np.diagonal(D)
        _Q = Q
    else:
        msg = "grid_matrix or D and Q unspecified."
        raise RuntimeError(msg)

    transformed_rots = np.zeros(rotations.shape, dtype='int_', order='C')
    is_compatible = ksp.snf_transform_rotations(
        transformed_rots,
        np.array(rotations, dtype='intc', order='C'),
        np.array(_D_diag, dtype='int_', order='C'),
        np.array(_Q, dtype='int_', order='C'))

    if is_compatible:
        return transformed_rots
    else:
        msg = "Grid generation matrix and rotation matrices are incompatible."
        raise RuntimeError(msg)


def get_all_grgrid_addresses(D_diag):
    """Return all grid addresses for mesh

    Parameters
    ----------
    D_diag : array_like
        Diagonal elements of D of Smith normal form.
        shape=(3,), dtype='int_'

    Returns
    -------
    grgrid_address : ndarray
        Genralized regular grid addresses of all grid points corresponding
        to D_diag.
        shape=(all_grid_points, 3), dtype='int_'

    """

    grgrid_address = np.zeros((np.prod(D_diag), 3), dtype='int_', order='C')
    ksp.all_grgrid_addresses(grgrid_address, np.array(D_diag, dtype='int_'))
    return grgrid_address


def get_double_grgrid_address(address, D_diag, PS=None):
    """Convert grid address plus shift to double-grid address

    address_double = 2 * address + shift, where values are reduced to
    be closet to 0, -mesh/2 < address_double <= mesh/2.

    Parameters
    ----------
    address : array_like
        Grid address.
        shape=(3,), dtype='int_'
    D_diag : array_like
        Diagonal elements of D of Smith normal form.
        shape=(3,), dtype='int_'
    PS : array_like, optional
        Half grid shifts after transformation by P of Smith normal form.
        Let half grid shifts along reciprocal basis vector directions be S,
        where s_i = 0 or 1, this array corresponds to np.dot(P, S).
        shape=(3,), dtype='int_'

    Returns
    -------
    address_double : ndarray
        Double-grid address.
        shape=(3,), dtype='intc'

    """

    address_double = np.zeros(3, dtype='int_', order='C')
    if PS is None:
        _PS = np.zeros(3, dtype='int_')
    else:
        _PS = np.array(PS, dtype='int_')
    ksp.double_grgrid_address(address_double,
                              np.array(address, dtype='int_'),
                              np.array(D_diag, dtype='int_'),
                              _PS)
    return address_double


def get_grgrid_index(address, D_diag):
    """Return grid point index of a single-grid address

    Parameters
    ----------
    address : array_like
        Single-grid address.
        shape=(3,), dtype='int_'
    D_diag : array_like
        Diagonal elements of D of Smith normal form.
        shape=(3,), dtype='int_'

    Returns
    -------
    grid_index : int
        Grid point index.

    """
    return ksp.grgrid_index(np.array(address, dtype='int_'),
                            np.array(D_diag, dtype='int_'))


def get_double_grgrid_index(address_double, D_diag, PS=None):
    """Return grid point index of a double-grid address

    Parameters
    ----------
    address_double : array_like
        Double-grid address.
        shape=(3,), dtype='int_'
    D_diag : array_like
        Diagonal elements of D of Smith normal form.
        shape=(3,), dtype='int_'
    PS : array_like, optional
        Half grid shifts after transformation by P of Smith normal form.
        Let half grid shifts along reciprocal basis vector directions be S,
        where s_i = 0 or 1, this array corresponds to np.dot(P, S).
        shape=(3,), dtype='int_'

    Returns
    -------
    grid_index : int
        Grid point index.

    """
    if PS is None:
        _PS = np.zeros(3, dtype='int_')
    else:
        _PS = np.array(PS, dtype='int_')
    return ksp.double_grgrid_index(np.array(address_double, dtype='int_'),
                                   np.array(D_diag, dtype='int_'),
                                   _PS)


def niggli_reduce(lattice, eps=1e-5):
    """Perform Niggli reduction

    Parameters
    ----------
    lattice : array_like
        Basis vectors in column vectors.
        shape=(3,3), dtype='double', order='C'
    eps : float, optional
        Tolerance parameter. Default is 1e-5.

    """

    red_lattice = np.zeros((3, 3), dtype='double', order='C')
    succeeded = ksp.niggli_reduce(red_lattice,
                                  np.array(lattice, dtype='double', order='C'),
                                  float(eps))
    if succeeded:
        return red_lattice
    else:
        msg = "Niggli reduction failed."
        raise RuntimeError(msg)
