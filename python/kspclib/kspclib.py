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
    """Return all grid addresses for mesh"""
    grid_address = np.zeros((np.prod(mesh), 3), dtype='intc', order='C')
    ksp.all_grid_addresses(grid_address, np.array(mesh, dtype='intc'))
    return grid_address


def get_grid_point_double_mesh(address_double, mesh):
    """Return grid point index of grid address of mesh"""
    return ksp.grid_point_double_mesh(np.array(address_double, dtype='intc'),
                                      np.array(mesh, dtype='intc'))


def get_grid_address_double_mesh(address, mesh, is_shift=None):
    """Return grid point index of grid address of mesh"""
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
