/* Copyright (C) 2020 Atsushi Togo */
/* All rights reserved. */

/* This file is part of kspclib. */

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
#include "kgengrid.h"
#include "mathfunc.h"
#include "snf3x3.h"

#define IDENTITY_TOL 1e-5

static void reduce_single_grid_address(long address[3], const long D_diag[3]);
static void reduce_double_grid_address(long address_double[3],
                                       const long D_diag[3]);
static size_t get_double_grid_point_index(const long address_double[3],
                                          const long D_diag[3],
                                          const long PS[3]);
static size_t get_grid_point_index(const long address[3],
                                   const long D_diag[3]);
static void get_all_grid_addresses(long grid_address[][3],
                                   const long D_diag[3]);

int kgg_get_snf3x3(long D_diag[3],
                   long P[3][3],
                   long Q[3][3],
                   MATCONST long A[3][3])
{
  int i, j, succeeded;
  long D[3][3];

  succeeded = 0;

  if (mat_get_determinant_l3(A) == 0) {
    goto err;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      D[i][j] = A[i][j];
    }
  }

  succeeded = snf3x3(D, P, Q);
  for (i = 0; i < 3; i++) {
    D_diag[i] = D[i][i];
  }

err:
  return succeeded;
}

/*----------------------------------------*/
/* Transform rotations by D(Q^-1)RQ(D^-1) */
/*----------------------------------------*/
/* transformed_rots : D(Q^-1)RQ(D^-1) */
/* rotations : [num_rot][3][3] */
/*    Defined as q' = Rq where q is in the reciprocal primitive basis */
/*    vectors. */
/* num_rot : Number of rotations */
int kgg_transform_rotations(long (*transformed_rots)[3][3],
                            MATCONST int (*rotations)[3][3],
                            const int num_rot,
                            const long D_diag[3],
                            MATCONST long Q[3][3])
{
  int i, j, k;
  double r[3][3], Q_double[3][3];

  /* Compute D(Q^-1)RQ(D^-1) by three steps */
  /* It is assumed that |det(Q)|=1 and Q^-1 has relatively small round-off */
  /* error, and we want to divide by D carefully. */
  /* 1. Compute (Q^-1)RQ */
  /* 2. Compute D(Q^-1)RQ */
  /* 3. Compute D(Q^-1)RQ(D^-1) */
  mat_cast_matrix_3l_to_3d(Q_double, Q);
  for (i = 0; i < num_rot; i++) {
    mat_get_similar_matrix_id3(r, rotations[i], Q_double, 0);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        r[j][k] *= D_diag[j];
        r[j][k] /= D_diag[k];
      }
    }
    mat_cast_matrix_3d_to_3l(transformed_rots[i], r);
    if (!mat_check_identity_matrix_ld3(transformed_rots[i], r, IDENTITY_TOL)) {
      return 0;
    }
  }

  return 1;
}

/* -------------------------------*/
/* Get all address in single grid */
/* -------------------------------*/
/* address : Single grid address. */
/* D_diag : Diagnal elements of D. */
void kgg_get_all_grid_addresses(long grid_address[][3], const long D_diag[3])
{
  get_all_grid_addresses(grid_address, D_diag);
}

/* -------------------------------------------------------*/
/* Get address in double grid from address in single grid */
/* -------------------------------------------------------*/
/* address_double : Double grid address. */
/* address : Single grid address. */
/* D_diag : Diagnal elements of D. */
/* PS : Shifts transformed by P. s_i is 0 or 1. */
void kgg_get_double_grid_address(long address_double[3],
                                 const long address[3],
                                 const long D_diag[3],
                                 const long PS[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    address_double[i] = address[i] * 2 + PS[i];
  }
  reduce_double_grid_address(address_double, D_diag);
}

/* -------------------------------------------------------*/
/* Get address in single grid from address in double grid */
/* -------------------------------------------------------*/
/* address : Single grid address. */
/* address_double : Double grid address. */
/* D_diag : Diagnal elements of D. */
/* PS : Shifts transformed by P. s_i is 0 or 1. */
void kgg_get_single_grid_address(long address[3],
                                 const long address_double[3],
                                 const long D_diag[3],
                                 const long PS[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    address[i] = (address_double[i] - PS[i]) / 2;
  }
  reduce_single_grid_address(address, D_diag);
}

/* -------------------------------------------------*/
/* Get grid point index from address in double grid */
/* -------------------------------------------------*/
/* address_double : Double grid address. */
/* D_diag : Diagnal elements of D. */
/* PS : Shifts transformed by P. s_i is 0 or 1. */
size_t kgg_get_double_grid_point(const long address_double[3],
                                 const long D_diag[3],
                                 const long PS[3])
{
  return get_double_grid_point_index(address_double, D_diag, PS);
}

static void reduce_single_grid_address(long address[3], const long D_diag[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    address[i] = mat_modulo_l(address[i], D_diag[i]);
  }
}

static void reduce_double_grid_address(long address_double[3],
                                       const long D_diag[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    address_double[i] = mat_modulo_l(address_double[i], 2 * D_diag[i]);
  }
}

static size_t get_double_grid_point_index(const long address_double[3],
                                          const long D_diag[3],
                                          const long PS[3])
{
  long address[3];

  kgg_get_single_grid_address(address,
                              address_double,
                              D_diag,
                              PS);
  return get_grid_point_index(address, D_diag);
}

/* See kgrid.h about GRID_ORDER_XYZ information. */
static size_t get_grid_point_index(const long address[3],
                                   const long D_diag[3])
{
#ifndef GRID_ORDER_XYZ
  return (address[2] * D_diag[0] * (size_t)(D_diag[1])
          + address[1] * D_diag[0] + address[0]);
#else
  return (address[0] * D_diag[1] * (size_t)(D_diag[2])
          + address[1] * D_diag[2] + address[2]);
#endif
}

static void get_all_grid_addresses(long grid_address[][3],
                                   const long D_diag[3])
{
  size_t i, j, k;
  size_t grid_point;
  long address[3];

  for (i = 0; i < (size_t)D_diag[0]; i++) {
    address[0] = i;
    for (j = 0; j < (size_t)D_diag[1]; j++) {
      address[1] = j;
      for (k = 0; k < (size_t)D_diag[2]; k++) {
        address[2] = k;
        grid_point = get_grid_point_index(address, D_diag);
        mat_copy_vector_l3(grid_address[grid_point], address);
      }
    }
  }
}
