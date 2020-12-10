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
#include "kgrid.h"
#include "kgengrid.h"
#include "mathfunc.h"
#include "snf3x3.h"
#include "tetrahedron_method.h"
#include "version.h"

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
  kgd_get_all_grid_addresses(grid_address, mesh);
}

size_t ksp_get_grid_point_double_mesh(const int address_double[3],
                                      const int mesh[3])
{
  return kgd_get_grid_point_double_mesh(address_double, mesh);
}

void ksp_get_grid_address_double_mesh(int address_double[3],
                                      const int address[3],
                                      const int mesh[3],
                                      const int is_shift[3])
{
  kgd_get_grid_address_double_mesh(address_double,
                                   address,
                                   mesh,
                                   is_shift);
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
  return kgg_get_snf3x3(D_diag, P, Q, A);
}

int ksp_snf_transform_rotations(long (*transformed_rots)[3][3],
                                KSPCONST int (*rotations)[3][3],
                                const int num_rot,
                                const long D_diag[3],
                                KSPCONST long Q[3][3])
{
  int succeeded;

  succeeded = kgg_transform_rotations(transformed_rots,
                                      rotations, num_rot, D_diag, Q);

  return succeeded;
}
