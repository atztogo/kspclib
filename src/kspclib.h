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

#include <stddef.h>
#include "tetrahedron_method.h"

int ksp_get_major_version(void);
int ksp_get_minor_version(void);
int ksp_get_micro_version(void);
void ksp_get_all_grid_addresses(int grid_address[][3], const int mesh[3]);
size_t ksp_get_grid_point_double_mesh(const int address_double[3],
                                      const int mesh[3]);
void ksp_get_grid_address_double_mesh(int address_double[3],
                                      const int address[3],
                                      const int mesh[3],
                                      const int is_shift[3]);
/* rec_lattice : column vectors */
void ksp_get_thm_relative_grid_addresses(int relative_grid_addresses[24][4][3],
                                         THMCONST double rec_lattice[3][3]);

#endif
