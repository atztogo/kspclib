/* Copyright (C) 2015 Atsushi Togo */
/* All rights reserved. */

/* This file was originally part of spglib and is part of kspclib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the kspclib project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "mathfunc.h"
#include "rgrid.h"

static void get_all_grid_addresses(int grid_address[][3], const int mesh[3]);
static long get_double_grid_index(const int address_double[3],
                                  const int mesh[3]);
static long get_grid_index_single_mesh(const int address[3],
                                       const int mesh[3]);
static void reduce_grid_address(int address[3], const int mesh[3]);
static void reduce_double_grid_address(int address[3], const int mesh[3]);

void rgd_get_all_grid_addresses(int grid_address[][3], const int mesh[3])
{
  get_all_grid_addresses(grid_address, mesh);
}

long rgd_get_double_grid_index(const int address_double[3],
                               const int mesh[3])
{
  return get_double_grid_index(address_double, mesh);
}

void rgd_get_double_grid_address(int address_double[3],
                                 const int address[3],
                                 const int mesh[3],
                                 const int is_shift[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    address_double[i] = address[i] * 2 + (is_shift[i] != 0);
  }
  reduce_double_grid_address(address_double, mesh);
}

static void get_all_grid_addresses(int grid_address[][3], const int mesh[3])
{
  int i, j, k;
  long grid_index;
  int address[3];

  for (i = 0; i < mesh[0]; i++) {
    address[0] = i;
    for (j = 0; j < mesh[1]; j++) {
      address[1] = j;
      for (k = 0; k < mesh[2]; k++) {
        address[2] = k;
        grid_index = get_grid_index_single_mesh(address, mesh);

        assert(mesh[0] * mesh[1] * mesh[2] > grid_index);

        grid_address[grid_index][0] = address[0];
        grid_address[grid_index][1] = address[1];
        grid_address[grid_index][2] = address[2];
        reduce_grid_address(grid_address[grid_index], mesh);
      }
    }
  }
}

static long get_double_grid_index(const int address_double[3],
                                  const int mesh[3])
{
  int i;
  int address[3];

  for (i = 0; i < 3; i++) {
    if (address_double[i] % 2 == 0) {
      address[i] = address_double[i] / 2;
    } else {
      address[i] = (address_double[i] - 1) / 2;
    }
    address[i] = mat_modulo_i(address[i], mesh[i]);
  }

  return get_grid_index_single_mesh(address, mesh);
}

static long get_grid_index_single_mesh(const int address[3],
                                       const int mesh[3])
{
#ifndef GRID_ORDER_XYZ
  return (address[2] * mesh[0] * (long)(mesh[1])
          + address[1] * mesh[0] + address[0]);
#else
  return (address[0] * mesh[1] * (long)(mesh[2])
          + address[1] * mesh[2] + address[2]);
#endif
}

static void reduce_grid_address(int address[3], const int mesh[3])
{
  int i;

  for (i = 0; i < 3; i++) {
#ifndef GRID_BOUNDARY_AS_NEGATIVE
    address[i] -= mesh[i] * (address[i] > mesh[i] / 2);
#else
    address[i] -= mesh[i] * (address[i] > (mesh[i] - 1) / 2);
#endif
  }
}

static void reduce_double_grid_address(int address[3], const int mesh[3])
{
  int i;

  for (i = 0; i < 3; i++) {
#ifndef GRID_BOUNDARY_AS_NEGATIVE
    address[i] -= 2 * mesh[i] * (address[i] > mesh[i]);
#else
    address[i] -= 2 * mesh[i] * (address[i] > mesh[i] - 1);
#endif
  }
}
