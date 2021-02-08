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

#ifndef __kspclib_H__
#define __kspclib_H__

#ifdef __cplusplus
extern "C" {
#endif

/* KSPCONST is used instead of 'const' so to avoid gcc warning. */
/* However there should be better way than this way.... */
#ifndef KSPCONST
#define KSPCONST
#endif

#include <stddef.h>

  int ksp_get_major_version(void);
  int ksp_get_minor_version(void);
  int ksp_get_micro_version(void);
  void ksp_get_all_grid_addresses(int grid_address[][3], const int mesh[3]);
  void ksp_get_double_grid_address(int address_double[3],
                                   const int address[3],
                                   const int mesh[3],
                                   const int is_shift[3]);
  long ksp_get_double_grid_index(const int address_double[3],
                                 const int mesh[3]);
/* rec_lattice : column vectors */
  void ksp_get_thm_relative_grid_addresses(int relative_grid_addresses[24][4][3],
                                           KSPCONST double rec_lattice[3][3]);
  double ksp_get_thm_integration_weight(const double omega,
                                        KSPCONST double tetrahedra_omegas[24][4],
                                        const char function);
  int ksp_get_snf3x3(long D_diag[3],
                     long P[3][3],
                     long Q[3][3],
                     KSPCONST long A[3][3]);
  int ksp_snf_transform_rotations(long (*transformed_rots)[3][3],
                                  KSPCONST long (*rotations)[3][3],
                                  const int num_rot,
                                  const long D_diag[3],
                                  KSPCONST long Q[3][3]);
  void ksp_get_all_grgrid_addresses(long grid_address[][3],
                                    const long D_diag[3]);
  void ksp_get_double_grgrid_address(long address_double[3],
                                     const long address[3],
                                     const long D_diag[3],
                                     const long PS[3]);
  long ksp_get_grgrid_index(const long address[3], const long D_diag[3]);
  long ksp_get_double_grgrid_index(const long address_double[3],
                                   const long D_diag[3],
                                   const long PS[3]);
  void ksp_get_grgrid_address_from_index(long address[3],
                                         const long grid_index,
                                         const long D_diag[3]);
  long ksp_rotate_grgrid_index(const long grid_index,
                               KSPCONST long rotation[3][3],
                               const long D_diag[3],
                               const long PS[3]);
  void ksp_get_ir_grgrid_map(long ir_grid_indices[],
                             KSPCONST long (*rotations)[3][3],
                             const int num_rot,
                             const long D_diag[3],
                             const long PS[3]);
  int ksp_niggli_reduce(double red_lattice[3][3],
                        KSPCONST double lattice[3][3],
                        const double eps);
  int ksp_get_reciprocal_point_group(long rec_rotations[48][3][3],
                                     KSPCONST long (*rotations)[3][3],
                                     const int num_rot,
                                     const int is_time_reversal);

#ifdef __cplusplus
}
#endif

#endif
