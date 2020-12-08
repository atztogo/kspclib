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

int kgg_get_snf3x3(long D[3][3],
                   long P[3][3],
                   long Q[3][3],
                   MATCONST long A[3][3])
{
  int i, j, succeeded;

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

err:
  return succeeded;
}

/* rotations->mat: (long*)[3][3] */
/*    Defined as q' = Rq */
int kgg_sanity_check_rotations(MATCONST MatINT *rotations,
                               MATCONST long D[3][3],
                               MATCONST long Q[3][3])
{
  int i, j, k;
  double r[3][3], Q_double[3][3];
  long r_long[3][3];

  /* Compute D(Q^-1)RQ(D^-1) by three steps */
  /* It is assumed that |det(Q)|=1 and Q^-1 has relatively small round-off */
  /* error, and we want to divide by D carefully. */
  /* 1. Compute (Q^-1)RQ */
  /* 2. Compute D(Q^-1)RQ */
  /* 3. Compute D(Q^-1)RQ(D^-1) */
  mat_cast_matrix_3l_to_3d(Q_double, Q);
  for (i = 0; i < rotations->size; i++) {
    mat_get_similar_matrix_id3(r, rotations->mat[i], Q_double, 0);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        r[j][k] *= D[j][j];
        r[j][k] /= D[k][k];
      }
    }
    mat_cast_matrix_3d_to_3l(r_long, r);
    if (! mat_check_identity_matrix_ld3(r_long, r, IDENTITY_TOL)) {
      return 0;
    }
  }

  return 1;
}
