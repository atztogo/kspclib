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


int kgg_get_snf3x3(SNF3x3 * snf, MATCONST long A[3][3])
{
  int i, j, succeeded;

  succeeded = 0;

  if (mat_get_determinant_l3(A) == 0) {
    goto err;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      snf->D[i][j] = A[i][j];
    }
  }

  succeeded = snf3x3(snf->D, snf->P, snf->Q);

err:
  return succeeded;
}

/* rotations->mat: (long*)[3][3] */
/*    Defined as q' = Rq */
int kgg_sanity_check_rotations(MATCONST SNF3x3 *snf,
                               MATCONST MatINT *rotations,
                               const double symprec)
{
  int i, j;
  double QD_inv[3][3], D_inv[3][3], r[3][3];
  long r_long[3][3];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i == j) {
        D_inv[i][i] = 1.0 / ((double) snf->D[i][i]);
      } else {
        D_inv[i][j] = 0;
      }
    }
  }

  /* D(Q^-1)RQ(D^-1) */
  mat_multiply_matrix_ld3(QD_inv, snf->Q, D_inv);
  for (i = 0; i < rotations->size; i++) {
    mat_get_similar_matrix_id3(r, rotations->mat[i], QD_inv, 0);
    mat_cast_matrix_3d_to_3l(r_long, r);
    if (! mat_check_identity_matrix_ld3(r_long, r, symprec)) {
      return 0;
    }
  }

  return 1;
}
