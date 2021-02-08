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

#include "kspclib.h"
#include "grgrid.h"
#include "mathfunc.h"
#include "niggli.h"
#include "rgrid.h"
#include "snf3x3.h"
#include "tetrahedron_method.h"
#include "version.h"

static int get_reciprocal_point_group(long rec_rotations[48][3][3],
                                      KSPCONST long (*rotations)[3][3],
                                      const int num_rot,
                                      const int is_time_reversal);

int ksp_get_major_version(void)
{
  return KSPCLIB_MAJOR_VERSION;
}

int ksp_get_minor_version(void)
{
  return KSPCLIB_MINOR_VERSION;
}

int ksp_get_micro_version(void)
{
  return KSPCLIB_MICRO_VERSION;
}

void ksp_get_all_grid_addresses(int grid_address[][3], const int mesh[3])
{
  rgd_get_all_grid_addresses(grid_address, mesh);
}

void ksp_get_double_grid_address(int address_double[3],
                                 const int address[3],
                                 const int mesh[3],
                                 const int is_shift[3])
{
  rgd_get_double_grid_address(address_double,
                              address,
                              mesh,
                              is_shift);
}

long ksp_get_double_grid_index(const int address_double[3],
                               const int mesh[3])
{
  return rgd_get_double_grid_index(address_double, mesh);
}

void ksp_get_thm_relative_grid_addresses(int relative_grid_addresses[24][4][3],
                                         KSPCONST double rec_lattice[3][3])
{
  thm_get_relative_grid_address(relative_grid_addresses, rec_lattice);
}

double ksp_get_thm_integration_weight(const double omega,
                                      KSPCONST double tetrahedra_omegas[24][4],
                                      const char function)
{
  return thm_get_integration_weight(omega, tetrahedra_omegas, function);
}

int ksp_get_snf3x3(long D_diag[3],
                   long P[3][3],
                   long Q[3][3],
                   KSPCONST long A[3][3])
{
  return grg_get_snf3x3(D_diag, P, Q, A);
}

int ksp_snf_transform_rotations(long (*transformed_rots)[3][3],
                                KSPCONST long (*rotations)[3][3],
                                const int num_rot,
                                const long D_diag[3],
                                KSPCONST long Q[3][3])
{
  int succeeded;

  succeeded = grg_transform_rotations(transformed_rots,
                                      rotations, num_rot, D_diag, Q);

  return succeeded;
}

void ksp_get_all_grgrid_addresses(long grid_address[][3],
                                  const long D_diag[3])
{
  grg_get_all_grid_addresses(grid_address, D_diag);
}

void ksp_get_double_grgrid_address(long address_double[3],
                                   const long address[3],
                                   const long D_diag[3],
                                   const long PS[3])
{
  grg_get_double_grid_address(address_double,
                              address,
                              D_diag,
                              PS);
}

long ksp_get_grgrid_index(const long address[3], const long D_diag[3])
{
  return grg_get_grid_index(address, D_diag);
}

long ksp_get_double_grgrid_index(const long address_double[3],
                                 const long D_diag[3],
                                 const long PS[3])
{
  return grg_get_double_grid_index(address_double,
                                   D_diag,
                                   PS);
}

void ksp_get_grgrid_address_from_index(long address[3],
                                       const long grid_index,
                                       const long D_diag[3])
{
  grg_get_grid_address_from_index(address, grid_index, D_diag);
}

long ksp_rotate_grgrid_index(const long grid_index,
                             KSPCONST long rotation[3][3],
                             const long D_diag[3],
                             const long PS[3])
{
  return grg_rotate_grid_index(grid_index, rotation, D_diag, PS);
}

void ksp_get_ir_grgrid_map(long ir_grid_indices[],
                           KSPCONST long (*rotations)[3][3],
                           const int num_rot,
                           const long D_diag[3],
                           const long PS[3])
{
  grg_get_ir_grid_map(ir_grid_indices,
                      rotations,
                      num_rot,
                      D_diag,
                      PS);
}

/* red_lattice, lattice : column vectors */
int ksp_niggli_reduce(double red_lattice[3][3],
                      KSPCONST double lattice[3][3],
                      const double eps)
{
  int i, j, succeeded;
  double lat[9];

  succeeded = 0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lat[i * 3 + j] = lattice[i][j];
    }
  }

  succeeded = niggli_reduce(lat, eps);

  if (succeeded) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        red_lattice[i][j] = lat[i * 3 + j];
      }
    }
  }

  return succeeded;
}

int ksp_get_reciprocal_point_group(long rec_rotations[48][3][3],
                                   KSPCONST long (*rotations)[3][3],
                                   const int num_rot,
                                   const int is_time_reversal)
{
  return get_reciprocal_point_group(rec_rotations,
                                    rotations,
                                    num_rot,
                                    is_time_reversal);
}

/* Extract unique rotation matrices and transpose them. */
/* When is_time_reversal == 1, inverse of the extracted matrices are */
/* included. */
/* Return 0 if failed */
static int get_reciprocal_point_group(long rec_rotations[48][3][3],
                                      KSPCONST long (*rotations)[3][3],
                                      const int num_rot,
                                      const int is_time_reversal)
{
  int i, j, num_rot_ret, inv_exist;
  KSPCONST long inversion[3][3] = {
    {-1, 0, 0 },
    { 0,-1, 0 },
    { 0, 0,-1 }
  };

  num_rot_ret = 0;
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < num_rot_ret; j++) {
      if (mat_check_identity_matrix_l3(rotations[i], rec_rotations[j])) {
        goto escape;
      }
    }
    if (num_rot_ret == 48) {
      goto err;
    }
    mat_copy_matrix_l3(rec_rotations[num_rot_ret], rotations[i]);
    num_rot_ret++;
  escape:
    ;
  }

  inv_exist = 0;
  if (is_time_reversal) {
    for (i = 0; i < num_rot_ret; i++) {
      if (mat_check_identity_matrix_l3(inversion, rec_rotations[i])) {
        inv_exist = 1;
        break;
      }
    }

    if (!inv_exist) {
      if (num_rot_ret > 24) {
        goto err;
      }

      for (i = 0; i < num_rot_ret; i++) {
        mat_multiply_matrix_l3(rec_rotations[num_rot_ret + i],
                               inversion, rec_rotations[i]);
      }
      num_rot_ret *= 2;
    }
  }

  for (i = 0; i < num_rot_ret; i++) {
    mat_transpose_matrix_l3(rec_rotations[i], rec_rotations[i]);
  }

  return num_rot_ret;
err:
  return 0;
}
