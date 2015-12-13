/* Copyright (C) 2014 Atsushi Togo */
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

/* * Neither the name of the phonopy project nor the names of its */
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
/* tetrahedron_method.c */
/* Copyright (C) 2014 Atsushi Togo */

#include "tetrahedron_method.h"
#include "kgrid.h"
#include <math.h>
#include <stdio.h>

#ifdef THMWARNING
#include <stdio.h>
#define warning_print(...) fprintf(stderr,__VA_ARGS__)
#else
#define warning_print(...)
#endif

/*      6-------7             */
/*     /|      /|             */
/*    / |     / |             */
/*   4-------5  |             */
/*   |  2----|--3             */
/*   | /     | /              */
/*   |/      |/	              */
/*   0-------1	              */
/*  		              */
/*  i: vec        neighbours  */
/*  0: O          1, 2, 4     */
/*  1: a          0, 3, 5     */
/*  2: b          0, 3, 6     */
/*  3: a + b      1, 2, 7     */
/*  4: c          0, 5, 6     */
/*  5: c + a      1, 4, 7     */
/*  6: c + b      2, 4, 7     */
/*  7: c + a + b  3, 5, 6     */


static int main_diagonals[4][3] = {{ 1, 1, 1},  /* 0-7 */
				   {-1, 1, 1},  /* 1-6 */
				   { 1,-1, 1},  /* 2-5 */
				   { 1, 1,-1}}; /* 3-4 */

static int db_relative_grid_address[4][24][4][3] = {
  {
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  1,  0}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  0,  1}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 1,  1,  0}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  1,  1}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 0,  0,  1}, { 1,  0,  1}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 0,  0,  1}, { 0,  1,  1}, { 1,  1,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  1,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0,  0,  1}, { 0,  1,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  0,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  0,  1}, { 1,  0,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  0,  1}, {-1, -1,  0}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  0,  1}, {-1, -1,  0}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  1,  0}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 1,  1,  0}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  0, -1}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  0, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0, -1, -1}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0, -1, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1, -1, -1}, { 0, -1, -1}, { 0,  0, -1}, },
    { { 0,  0,  0}, {-1, -1, -1}, { 0, -1, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1, -1, -1}, {-1,  0, -1}, { 0,  0, -1}, },
    { { 0,  0,  0}, {-1, -1, -1}, {-1,  0, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1, -1, -1}, {-1, -1,  0}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1, -1, -1}, {-1, -1,  0}, {-1,  0,  0}, },
  },
  {
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  1,  0}, { 0,  1,  1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  0,  1}, { 0,  1,  1}, },
    { { 0,  0,  0}, {-1,  1,  0}, {-1,  1,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  0,  1}, {-1,  1,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  1,  0}, {-1,  1,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  1,  1}, { 0,  1,  1}, },
    { { 0,  0,  0}, {-1,  0,  1}, { 0,  0,  1}, {-1,  1,  1}, },
    { { 0,  0,  0}, { 0,  0,  1}, {-1,  1,  1}, { 0,  1,  1}, },
    { { 0,  0,  0}, { 0,  0,  1}, { 0, -1,  0}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  0,  1}, { 1, -1,  0}, },
    { { 0,  0,  0}, {-1,  0,  1}, { 0, -1,  0}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  0,  1}, { 0,  0,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  0, -1}, { 1,  0, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  1,  0}, { 1,  0, -1}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  0, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  1,  0}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 0, -1, -1}, { 1, -1, -1}, { 0,  0, -1}, },
    { { 0,  0,  0}, { 0, -1, -1}, { 1, -1, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 1, -1, -1}, { 0,  0, -1}, { 1,  0, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1, -1, -1}, { 1,  0, -1}, },
    { { 0,  0,  0}, { 1, -1, -1}, { 0, -1,  0}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1, -1, -1}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 0, -1, -1}, { 0,  0, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0, -1, -1}, { 0, -1,  0}, {-1,  0,  0}, },
  },
  {
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  1,  0}, { 1,  0,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1,  0,  1}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  0,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  1,  0}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 1, -1,  1}, { 0, -1,  0}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 0, -1,  1}, { 1, -1,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1, -1,  1}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1, -1,  1}, { 1,  0,  1}, },
    { { 0,  0,  0}, { 0, -1,  1}, { 1, -1,  1}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 1, -1,  1}, { 0,  0,  1}, { 1,  0,  1}, },
    { { 0,  0,  0}, { 0, -1,  1}, { 0, -1,  0}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0, -1,  1}, { 0,  0,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  0, -1}, { 0,  1, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  1,  0}, { 0,  1, -1}, },
    { { 0,  0,  0}, {-1,  0, -1}, { 0,  0, -1}, {-1,  1, -1}, },
    { { 0,  0,  0}, {-1,  0, -1}, {-1,  1, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0,  0, -1}, {-1,  1, -1}, { 0,  1, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  1, -1}, { 0,  1, -1}, },
    { { 0,  0,  0}, {-1,  1,  0}, {-1,  1, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1,  1,  0}, { 0,  1,  0}, {-1,  1, -1}, },
    { { 0,  0,  0}, { 0,  0, -1}, { 0, -1,  0}, { 1, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0,  0, -1}, { 1, -1,  0}, },
    { { 0,  0,  0}, {-1,  0, -1}, { 0,  0, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1,  0, -1}, { 0, -1,  0}, {-1,  0,  0}, },
  },
  {
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  1,  0}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 1,  1,  0}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  0,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0,  1,  0}, {-1,  0,  1}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0, -1,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 0, -1,  1}, { 0,  0,  1}, },
    { { 0,  0,  0}, {-1, -1,  1}, {-1, -1,  0}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1, -1,  1}, {-1, -1,  0}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1, -1,  1}, { 0, -1,  1}, { 0, -1,  0}, },
    { { 0,  0,  0}, {-1, -1,  1}, {-1,  0,  1}, {-1,  0,  0}, },
    { { 0,  0,  0}, {-1, -1,  1}, { 0, -1,  1}, { 0,  0,  1}, },
    { { 0,  0,  0}, {-1, -1,  1}, {-1,  0,  1}, { 0,  0,  1}, },
    { { 0,  0,  0}, { 0,  0, -1}, { 1,  0, -1}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 0,  0, -1}, { 0,  1, -1}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  0, -1}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  1, -1}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  1,  0}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 1,  1,  0}, { 1,  1, -1}, },
    { { 0,  0,  0}, { 0,  0, -1}, { 0,  1, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0,  1,  0}, { 0,  1, -1}, {-1,  0,  0}, },
    { { 0,  0,  0}, { 0,  0, -1}, { 1,  0, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 1,  0,  0}, { 1,  0, -1}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  0, -1}, {-1, -1,  0}, { 0, -1,  0}, },
    { { 0,  0,  0}, { 0,  0, -1}, {-1, -1,  0}, {-1,  0,  0}, },
  },
};

static void get_integration_weight_at_omegas(const int kpoint, double *integration_weights,
				 const int num_omegas,
				 const double *omegas,
				 THMCONST double tetrahedra_omegas[24][4],
				 const char function);
static double get_integration_weight(const int kpoint, const double omega,
				     THMCONST double tetrahedra_omegas[24][4],
				     const char function,
				     int bloechl);
static double get_vertex_integration_weight(const double omega,
					    const double v[4],
					    const int pos0,
					    double (*gn)(const int,
							 const double,
							 const double[4]),
					    double (*IJ)(const int,
							 const int,
							 const double,
							 const double[4]));
static int get_main_diagonal(THMCONST double rec_lattice[3][3]);
static int sort_omegas(double v[4]);
static double norm_squared_d3(const double a[3]);
static void multiply_matrix_vector_di3(double v[3],
				       THMCONST double a[3][3],
				       const int b[3]);
static double _f(const int n,
		 const int m,
		 const double omega,
		 const double vertices_omegas[4]);
static double _delta(const int n,
		     const int m,
		     const double vertices_omegas[4]);
static double _J(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4]);
static double _I(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4]);
static double _Ider(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4]);
static double _n(const int i,
		 const double omega,
		 const double vertices_omegas[4]);
static double _g(const int i,
		 const double omega,
		 const double vertices_omegas[4]);
static double _gder(const int i,
		 const double omega,
		 const double vertices_omegas[4]);
static double _n_0(void);
static double _n_1(const double omega,
		   const double vertices_omegas[4]);
static double _n_2(const double omega,
		   const double vertices_omegas[4]);
static double _n_3(const double omega,
		   const double vertices_omegas[4]);
static double _n_4(void);
static double _g_0(void);
static double _g_1(const double omega,
		   const double vertices_omegas[4]);
static double _g_2(const double omega,
		   const double vertices_omegas[4]);
static double _g_3(const double omega,
		   const double vertices_omegas[4]);
static double _g_4(void);
static double _gder_0(void);
static double _gder_1(const double omega,
		   const double vertices_omegas[4]);
static double _gder_2(const double omega,
		   const double vertices_omegas[4]);
static double _gder_3(const double omega,
		   const double vertices_omegas[4]);
static double _gder_4(void);
static double _J_0(void);
static double _J_10(const double omega,
		    const double vertices_omegas[4]);
static double _J_11(const double omega,
		    const double vertices_omegas[4]);
static double _J_12(const double omega,
		    const double vertices_omegas[4]);
static double _J_13(const double omega,
		    const double vertices_omegas[4]);
static double _J_20(const double omega,
		    const double vertices_omegas[4]);
static double _J_21(const double omega,
		    const double vertices_omegas[4]);
static double _J_22(const double omega,
		    const double vertices_omegas[4]);
static double _J_23(const double omega,
		    const double vertices_omegas[4]);
static double _J_30(const double omega,
		    const double vertices_omegas[4]);
static double _J_31(const double omega,
		    const double vertices_omegas[4]);
static double _J_32(const double omega,
		    const double vertices_omegas[4]);
static double _J_33(const double omega,
		    const double vertices_omegas[4]);
static double _J_4(void);
static double _I_0(void);
static double _I_10(const double omega,
		    const double vertices_omegas[4]);
static double _I_11(const double omega,
		    const double vertices_omegas[4]);
static double _I_12(const double omega,
		    const double vertices_omegas[4]);
static double _I_13(const double omega,
		    const double vertices_omegas[4]);
static double _I_20(const double omega,
		    const double vertices_omegas[4]);
static double _I_21(const double omega,
		    const double vertices_omegas[4]);
static double _I_22(const double omega,
		    const double vertices_omegas[4]);
static double _I_23(const double omega,
		    const double vertices_omegas[4]);
static double _I_30(const double omega,
		    const double vertices_omegas[4]);
static double _I_31(const double omega,
		    const double vertices_omegas[4]);
static double _I_32(const double omega,
		    const double vertices_omegas[4]);
static double _I_33(const double omega,
		    const double vertices_omegas[4]);
static double _I_4(void);

static double _Ider_0(void);
static double _Ider_10(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_11(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_12(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_13(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_20(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_21(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_22(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_23(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_30(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_31(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_32(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_33(const double omega,
		    const double vertices_omegas[4]);
static double _Ider_4(void);


void thm_get_relative_grid_address(int relative_grid_address[24][4][3],
				   THMCONST double rec_lattice[3][3])
{
  int i, j, k, main_diag_index;

  main_diag_index = get_main_diagonal(rec_lattice);
 
  for (i = 0; i < 24; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 3; k++) {
	relative_grid_address[i][j][k] =
	  db_relative_grid_address[main_diag_index][i][j][k];
      }
    }
  }
}

void thm_get_all_relative_grid_address(int relative_grid_address[4][24][4][3])
{
  int i, j, k, main_diag_index;
  
  for (main_diag_index = 0; main_diag_index < 4; main_diag_index++) {
    for (i = 0; i < 24; i++) {
      for (j = 0; j < 4; j++) {
	for (k = 0; k < 3; k++) {
	  relative_grid_address[main_diag_index][i][j][k] =
	    db_relative_grid_address[main_diag_index][i][j][k];
	}
      }
    }
  }
}

double thm_get_integration_weight(const int kpoint, const double omega,
				  THMCONST double tetrahedra_omegas[24][4],
				  const char function, int bloechl)
{
  return get_integration_weight(kpoint, omega,
				tetrahedra_omegas,
				function,
				bloechl);
}

void
thm_get_integration_weight_at_omegas(const int kpoint, double *integration_weights,
				     const int num_omegas,
				     const double *omegas,
				     THMCONST double tetrahedra_omegas[24][4],
				     const char function)
{
  get_integration_weight_at_omegas(kpoint, integration_weights,
				   num_omegas,
				   omegas,
				   tetrahedra_omegas,
				   function);
}

void thm_get_neighboring_grid_points(int neighboring_grid_points[],
				     const int grid_point,
				     THMCONST int relative_grid_address[][3],
				     const int num_relative_grid_address,
				     const int mesh[3],
				     THMCONST int bz_grid_address[][3],
				     const int bz_map[])
{
  int bzmesh[3], address_double[3], bz_address_double[3];
  int i, j, bz_gp;

  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2;
  }
  for (i = 0; i < num_relative_grid_address; i++) {
    for (j = 0; j < 3; j++) {
      address_double[j] = (bz_grid_address[grid_point][j] +
			   relative_grid_address[i][j]) * 2;
      bz_address_double[j] = address_double[j];
    }
    bz_gp = bz_map[kgd_get_grid_point_double_mesh(bz_address_double, bzmesh)];
    if (bz_gp == -1) {
      neighboring_grid_points[i] =
	kgd_get_grid_point_double_mesh(address_double, mesh);
    } else {
      neighboring_grid_points[i] = bz_gp;
    }
  }
}

static void
get_integration_weight_at_omegas(const int kpoint, double *integration_weights,
				 const int num_omegas,
				 const double *omegas,
				 THMCONST double tetrahedra_omegas[24][4],
				 const char function)
{
  int i;

#pragma omp parallel for
  for (i = 0; i < num_omegas; i++) {
    integration_weights[i] = get_integration_weight(kpoint,omegas[i],
						    tetrahedra_omegas,
						    function,
						    0);
  }
}

static double get_integration_weight(const int kpoint, const double omega,
				     THMCONST double tetrahedra_omegas[24][4],
				     const char function,
				     int bloechl)
{
  int i, j, pos0;
  double sum;
  double v[4];
  double (*gn)(const int, const double, const double[4]);
  double (*IJ)(const int, const int, const double, const double[4]);

  double tmpval;
  double standard, correction;
  
  if (function == 'I') {
    gn = _g;
    IJ = _I;
  } else {
    gn = _n;
    IJ = _J;
  }

  sum = 0;
  for (i = 0; i < 24; i++) {
    for (j = 0; j < 4; j++) {
      v[j] = tetrahedra_omegas[i][j];
    }
    pos0 = sort_omegas(v);
    standard=get_vertex_integration_weight(omega, v, pos0, gn, IJ);
    sum += standard;
    if (bloechl) {
      if (function=='J') {
        correction = (1.0 / 40 *
		       get_vertex_integration_weight(omega, v, pos0, _g, _I) *
		       (v[0] + v[1] + v[2] + v[3] - v[pos0] * 4));
      }
      else {
	// call twice, product rule for the derivative
	// or no, pretty sure we do not need to do this...
	// the correction is added to J(omega), if we take the derivative
	// we just get I(omega)+der(correction), where der(correction)
	// only considers the derivative of the density of states
	// (otherwise \sum I probably does not sum to one for one
	// tetrahedron if F=1...)
	correction=(1.0 / 40 *
		    get_vertex_integration_weight(omega, v, pos0, _gder, _I) *
		    (v[0] + v[1] + v[2] + v[3] - v[pos0] * 4));
      }
      sum += correction;
    }
  }

  return sum / 6;
}

static double get_vertex_integration_weight(const double omega,
					    const double v[4],
					    const int pos0,
					    double (*gn)(const int,
							 const double,
							 const double[4]),
					    double (*IJ)(const int,
							 const int,
							 const double,
							 const double[4]))
{
  double sum;

  sum = 0;
  if (omega < v[0]) {
    sum += IJ(0, pos0, omega, v) * gn(0, omega, v);
  } else {
    if (v[0] < omega && omega < v[1]) {
      sum += IJ(1, pos0, omega, v) * gn(1, omega, v);
    } else {
      if (v[1] < omega && omega < v[2]) {
	sum += IJ(2, pos0, omega, v) * gn(2, omega, v);
      } else {
	if (v[2] < omega && omega < v[3]) {
	  sum += IJ(3, pos0, omega, v) * gn(3, omega, v);
	} else {
	  if (v[3] < omega) {
	    sum += IJ(4, pos0, omega, v) * gn(4, omega, v);
	  }
	}
      }
    }
  }

  return sum;
}

static int sort_omegas(double v[4])
{
  int i;
  double w[4];

  /* Variable i gives the position of the vertex focused */
  /* (index0 before sorting) after sorting the vertices by incresing energy. */

  i = 0;
  
  if (v[0] > v[1]) {
    w[0] = v[1];
    w[1] = v[0];
    i = 1;
  } else {
    w[0] = v[0];
    w[1] = v[1];
  }

  if (v[2] > v[3]) {
    w[2] = v[3];
    w[3] = v[2];
  } else {
    w[2] = v[2];
    w[3] = v[3];
  }

  if (w[0] > w[2]) {
    v[0] = w[2];
    v[1] = w[0];
    if (i == 0) {
      i = 4;
    }
  } else {
    v[0] = w[0];
    v[1] = w[2];
  }

  if (w[1] > w[3]) {
    v[3] = w[1];
    v[2] = w[3];
    if (i == 1) {
      i = 3;
    }
  } else {
    v[3] = w[3];
    v[2] = w[1];
    if (i == 1) {
      i = 5;
    }
  }

  if (v[1] > v[2]) {
    w[1] = v[1];
    v[1] = v[2];
    v[2] = w[1];
    if (i == 4) {
      i = 2;
    }
    if (i == 5) {
      i = 1;
    }
  } else {
    if (i == 4) {
      i = 1;
    }
    if (i == 5) {
      i = 2;
    }
  }
  return i;
}

static int get_main_diagonal(THMCONST double rec_lattice[3][3])
{
  int i, shortest;
  double length, min_length;
  double main_diag[3];

  shortest = 0;
  multiply_matrix_vector_di3(main_diag, rec_lattice, main_diagonals[0]);
  min_length = norm_squared_d3(main_diag);
  for (i = 1; i < 4; i++) {
    multiply_matrix_vector_di3(main_diag, rec_lattice, main_diagonals[i]);
    length = norm_squared_d3(main_diag);
    if (min_length > length) {
      min_length = length;
      shortest = i;
    }
  }
  return shortest;
}

static double norm_squared_d3(const double a[3])
{
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

static void multiply_matrix_vector_di3(double v[3],
				       THMCONST double a[3][3],
				       const int b[3])
{
  int i;
  double c[3];

  for (i = 0; i < 3; i++) {
    c[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
  }

  for (i = 0; i < 3; i++) {
    v[i] = c[i];
  }
}

static double _f(const int n,
		 const int m,
		 const double omega,
		 const double vertices_omegas[4])
{
  return ((omega - vertices_omegas[m]) /
	  (vertices_omegas[n] - vertices_omegas[m]));
}

static double _delta(const int n,
		     const int m,
		     const double vertices_omegas[4])
{
  return vertices_omegas[n]-vertices_omegas[m];
}

static double _J(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _J_0();
  case 1:
    switch (pos0) {
    case 0:
      return _J_10(omega, vertices_omegas);
    case 1:
      return _J_11(omega, vertices_omegas);
    case 2:
      return _J_12(omega, vertices_omegas);
    case 3:
      return _J_13(omega, vertices_omegas);
    }
  case 2:
    switch (pos0) {
    case 0:
      return _J_20(omega, vertices_omegas);
    case 1:
      return _J_21(omega, vertices_omegas);
    case 2:
      return _J_22(omega, vertices_omegas);
    case 3:
      return _J_23(omega, vertices_omegas);
    }
  case 3:
    switch (pos0) {
    case 0:
      return _J_30(omega, vertices_omegas);
    case 1:
      return _J_31(omega, vertices_omegas);
    case 2:
      return _J_32(omega, vertices_omegas);
    case 3:
      return _J_33(omega, vertices_omegas);
    }
  case 4:
    return _J_4();
  }

  warning_print("******* Warning *******\n");
  warning_print(" J is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}


static double _I(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _I_0();
  case 1:
    switch (pos0) {
    case 0:
      return _I_10(omega, vertices_omegas);
    case 1:
      return _I_11(omega, vertices_omegas);
    case 2:
      return _I_12(omega, vertices_omegas);
    case 3:
      return _I_13(omega, vertices_omegas);
    }
  case 2:
    switch (pos0) {
    case 0:
      return _I_20(omega, vertices_omegas);
    case 1:
      return _I_21(omega, vertices_omegas);
    case 2:
      return _I_22(omega, vertices_omegas);
    case 3:
      return _I_23(omega, vertices_omegas);
    }
  case 3:
    switch (pos0) {
    case 0:
      return _I_30(omega, vertices_omegas);
    case 1:
      return _I_31(omega, vertices_omegas);
    case 2:
      return _I_32(omega, vertices_omegas);
    case 3:
      return _I_33(omega, vertices_omegas);
    }
  case 4:
    return _I_4();
  }

  warning_print("******* Warning *******\n");
  warning_print(" I is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}

static double _Ider(const int i,
		 const int pos0,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _Ider_0();
  case 1:
    switch (pos0) {
    case 0:
      return _Ider_10(omega, vertices_omegas);
    case 1:
      return _Ider_11(omega, vertices_omegas);
    case 2:
      return _Ider_12(omega, vertices_omegas);
    case 3:
      return _Ider_13(omega, vertices_omegas);
    }
  case 2:
    switch (pos0) {
    case 0:
      return _Ider_20(omega, vertices_omegas);
    case 1:
      return _Ider_21(omega, vertices_omegas);
    case 2:
      return _Ider_22(omega, vertices_omegas);
    case 3:
      return _Ider_23(omega, vertices_omegas);
    }
  case 3:
    switch (pos0) {
    case 0:
      return _Ider_30(omega, vertices_omegas);
    case 1:
      return _Ider_31(omega, vertices_omegas);
    case 2:
      return _Ider_32(omega, vertices_omegas);
    case 3:
      return _Ider_33(omega, vertices_omegas);
    }
  case 4:
    return _Ider_4();
  }

  warning_print("******* Warning *******\n");
  warning_print(" Ider is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}

static double _n(const int i,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _n_0();
  case 1:
    return _n_1(omega, vertices_omegas);
  case 2:
    return _n_2(omega, vertices_omegas);
  case 3:
    return _n_3(omega, vertices_omegas);
  case 4:
    return _n_4();
  }
  
  warning_print("******* Warning *******\n");
  warning_print(" n is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}

static double _g(const int i,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _g_0();
  case 1:
    return _g_1(omega, vertices_omegas);
  case 2:
    return _g_2(omega, vertices_omegas);
  case 3:
    return _g_3(omega, vertices_omegas);
  case 4:
    return _g_4();
  }
  
  warning_print("******* Warning *******\n");
  warning_print(" g is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}

static double _gder(const int i,
		 const double omega,
		 const double vertices_omegas[4])
{
  switch (i) {
  case 0:
    return _gder_0();
  case 1:
    return _gder_1(omega, vertices_omegas);
  case 2:
    return _gder_2(omega, vertices_omegas);
  case 3:
    return _gder_3(omega, vertices_omegas);
  case 4:
    return _gder_4();
  }
  
  warning_print("******* Warning *******\n");
  warning_print(" gder is something wrong. \n");
  warning_print("******* Warning *******\n");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

  return 0;
}

/* omega < omega1 */
static double _n_0(void)
{
  return 0.0;
}

/* omega1 < omega < omega2 */
static double _n_1(const double omega,
		   const double vertices_omegas[4])
{
  return (_f(1, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(3, 0, omega, vertices_omegas));
}

/* omega2 < omega < omega3 */
static double _n_2(const double omega,
		   const double vertices_omegas[4])
{
  return (_f(3, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas));
}
            
/* omega2 < omega < omega3 */
static double _n_3(const double omega,
		   const double vertices_omegas[4])
{
  return (1.0 -
	  _f(0, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas));
}

/* omega4 < omega */
static double _n_4(void)
{
  return 1.0;
}

/* omega < omega1 */
static double _g_0(void)
{
  return 0.0;
}

/* omega1 < omega < omega2 */
static double _g_1(const double omega,
		   const double vertices_omegas[4])
{
  return (3 *
	  _f(1, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) /
	  _delta(3,0,vertices_omegas));
}

/* omega2 < omega < omega3 */
static double _g_2(const double omega,
		   const double vertices_omegas[4])
{
  return (3 /
	  _delta(3, 0, vertices_omegas) *
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas)));
}

/* omega3 < omega < omega4 */
static double _g_3(const double omega,
		   const double vertices_omegas[4])
{
  return (3 *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas) /
	  _delta(3, 0, vertices_omegas));
}

/* omega4 < omega */
static double _g_4(void)
{
  return 0.0;
}

/* omega < omega1 */
static double _gder_0(void)
{
  return 0.0;
}

/* omega1 < omega < omega2 */
static double _gder_1(const double omega,
		   const double vertices_omegas[4])
{
  return (3 / _delta(3, 0, vertices_omegas) *
	  (_f(2, 0, omega, vertices_omegas) /
	   _delta(1, 0, vertices_omegas) +
	   _f(1, 0, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas)));
}

/* omega2 < omega < omega3 */
static double _gder_2(const double omega,
		   const double vertices_omegas[4])
{
  return (3 / _delta(3, 0, vertices_omegas) *
	  (_f(2, 0, omega, vertices_omegas)/_delta(1, 2, vertices_omegas) +
	   _f(1, 2, omega, vertices_omegas)/_delta(2, 0, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas)/_delta(2, 1, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas)/_delta(1, 3, vertices_omegas)));
}

/* omega3 < omega < omega4 */
static double _gder_3(const double omega,
		   const double vertices_omegas[4])
{
  return (3 / _delta(3, 0, vertices_omegas) *
	  (_f(2, 3, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 3, vertices_omegas)));
}

/* omega4 < omega */
static double _gder_4(void)
{
  return 0.0;
}

static double _J_0(void)
{
  return 0.0;
}

static double _J_10(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0 +
	  _f(0, 1, omega, vertices_omegas) +
	  _f(0, 2, omega, vertices_omegas) +
	  _f(0, 3, omega, vertices_omegas)) / 4;
}

static double _J_11(const double omega,
		    const double vertices_omegas[4])
{
  return _f(1, 0, omega, vertices_omegas) / 4;
}

static double _J_12(const double omega,
		    const double vertices_omegas[4])
{
  return _f(2, 0, omega, vertices_omegas) / 4;
}

static double _J_13(const double omega,
		    const double vertices_omegas[4])
{
  return _f(3, 0, omega, vertices_omegas) / 4;
}

static double _J_20(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (1.0 +
	   _f(0, 3, omega, vertices_omegas)) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) *
	  (1.0 +
	   _f(0, 3, omega, vertices_omegas) +
	   _f(0, 2, omega, vertices_omegas))) / 4 / _n_2(omega, vertices_omegas);
}

static double _J_21(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (1.0 +
	   _f(1, 3, omega, vertices_omegas) +
	   _f(1, 2, omega, vertices_omegas)) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (_f(1, 3, omega, vertices_omegas) +
	   _f(1, 2, omega, vertices_omegas)) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas)) / 4 / _n_2(omega, vertices_omegas);
}

static double _J_22(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) *
	  (_f(2, 1, omega, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas))) / 4 / _n_2(omega, vertices_omegas);
}

static double _J_23(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 1, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  _f(3, 1, omega, vertices_omegas) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (_f(3, 1, omega, vertices_omegas) +
	   _f(3, 0, omega, vertices_omegas)) +
	  _f(3, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) *
	  _f(3, 0, omega, vertices_omegas)) / 4 / _n_2(omega, vertices_omegas);
}

static double _J_30(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0 -
	  _f(0, 3, omega, vertices_omegas) *
	  _f(0, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas)) / 4 / _n_3(omega, vertices_omegas);
}

static double _J_31(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0 -
	  _f(0, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas)) / 4 / _n_3(omega, vertices_omegas);
}

static double _J_32(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0 +
	  _f(0, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas)) / 4 / _n_3(omega, vertices_omegas);
}

static double _J_33(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0 -
	  _f(0, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 3, omega, vertices_omegas) *
	  (1.0 +
	   _f(3, 0, omega, vertices_omegas) +
	   _f(3, 1, omega, vertices_omegas) +
	   _f(3, 2, omega, vertices_omegas))) / 4 / _n_3(omega, vertices_omegas);
}

static double _J_4(void)
{
  return 0.25;
}

static double _I_0(void)
{
  return 0.0;
}

static double _I_10(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(0, 1, omega, vertices_omegas) +
	  _f(0, 2, omega, vertices_omegas) +
	  _f(0, 3, omega, vertices_omegas)) / 3;
}

static double _I_11(const double omega,
		    const double vertices_omegas[4])
{
  return _f(1, 0, omega, vertices_omegas) / 3;
}

static double _I_12(const double omega,
		    const double vertices_omegas[4])
{
  return _f(2, 0, omega, vertices_omegas) / 3;
}

static double _I_13(const double omega,
		    const double vertices_omegas[4])
{
  return _f(3, 0, omega, vertices_omegas) / 3;
}

static double _I_20(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(0, 3, omega, vertices_omegas) +
	  _f(0, 2, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _I_21(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(1, 2, omega, vertices_omegas) +
	  _f(1, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _I_22(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(2, 1, omega, vertices_omegas) +
	  _f(2, 0, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}
            
static double _I_23(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 0, omega, vertices_omegas) +
	  _f(3, 1, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _I_30(const double omega,
		    const double vertices_omegas[4])
{
  return _f(0, 3, omega, vertices_omegas) / 3;
}

static double _I_31(const double omega,
		    const double vertices_omegas[4])
{
  return _f(1, 3, omega, vertices_omegas) / 3;
}

static double _I_32(const double omega,
		    const double vertices_omegas[4])
{
  return _f(2, 3, omega, vertices_omegas) / 3;
}

static double _I_33(const double omega,
		    const double vertices_omegas[4])
{
  return (_f(3, 0, omega, vertices_omegas) +
	  _f(3, 1, omega, vertices_omegas) +
	  _f(3, 2, omega, vertices_omegas)) / 3;
}

static double _I_4(void)
{
  return 0.0;
}

static double _Ider_0(void)
{
  return 0.0;
}

static double _Ider_10(const double omega,
		    const double vertices_omegas[4])
{
  return ( 1.0 / _delta(0, 1, vertices_omegas) +
	   1.0 / _delta(0, 2, vertices_omegas) +
	   1.0 / _delta(0, 3, vertices_omegas)) / 3;
}

static double _Ider_11(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(1, 0, vertices_omegas) / 3;
}

static double _Ider_12(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(2, 0, vertices_omegas) / 3;
}

static double _Ider_13(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(3, 0, vertices_omegas) / 3;
}

static double _Ider_20(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0/_delta(0, 3, vertices_omegas) -
	  _f(0, 2, omega, vertices_omegas) *
	  _f(2, 0, omega, vertices_omegas) *
	  _f(1, 2, omega, vertices_omegas) *
	  (_f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  pow(_f(1, 2, omega, vertices_omegas) *
	      _f(2, 0, omega, vertices_omegas) +
	      _f(2, 1, omega, vertices_omegas) *
	      _f(1, 3, omega, vertices_omegas),2.0) +
	  (_f(2, 0, omega, vertices_omegas) *
	   _f(1, 2, omega, vertices_omegas) /
	   _delta(0, 2, vertices_omegas) +
	   _f(0, 2, omega, vertices_omegas) *
	   _f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(0, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas)) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _Ider_21(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0/_delta(1, 2, vertices_omegas) -
	  _f(1, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (_f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  pow(_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	      _f(1, 3, omega, vertices_omegas),2.0) +
	  (2.0*_f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _Ider_22(const double omega,
		    const double vertices_omegas[4])
{
  return (1.0/_delta(1, 2, vertices_omegas) -
	  _f(1, 3, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (_f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  pow(_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	      _f(1, 3, omega, vertices_omegas),2.0) +
	  (2.0*_f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas)) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}
            
static double _Ider_23(const double omega,
		    const double vertices_omegas[4])
{
    return (1.0/_delta(3, 0, vertices_omegas) -
	  _f(3, 1, omega, vertices_omegas) *
	  _f(1, 3, omega, vertices_omegas) *
	  _f(2, 1, omega, vertices_omegas) *
	  (_f(1, 2, omega, vertices_omegas) /
	   _delta(2, 0, vertices_omegas) +
	   _f(2, 0, omega, vertices_omegas) /
	   _delta(1, 2, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  pow(_f(1, 2, omega, vertices_omegas) *
	      _f(2, 0, omega, vertices_omegas) +
	      _f(2, 1, omega, vertices_omegas) *
	      _f(1, 3, omega, vertices_omegas),2.0) +
	  (_f(1, 3, omega, vertices_omegas) *
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(3, 1, vertices_omegas) +
	   _f(3, 1, omega, vertices_omegas) *
	   _f(2, 1, omega, vertices_omegas) /
	   _delta(1, 3, vertices_omegas) +
	   _f(3, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas) /
	   _delta(2, 1, vertices_omegas)) /
	  (_f(1, 2, omega, vertices_omegas) *
	   _f(2, 0, omega, vertices_omegas) +
	   _f(2, 1, omega, vertices_omegas) *
	   _f(1, 3, omega, vertices_omegas))) / 3;
}

static double _Ider_30(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(0, 3, vertices_omegas) / 3;
}

static double _Ider_31(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(1, 3, vertices_omegas) / 3;
}

static double _Ider_32(const double omega,
		    const double vertices_omegas[4])
{
  return 1.0 / _delta(2, 3, vertices_omegas) / 3;
}

static double _Ider_33(const double omega,
		    const double vertices_omegas[4])
{
  return ( 1.0 / _delta(3, 0, vertices_omegas) +
	   1.0 / _delta(3, 1, vertices_omegas) +
	   1.0 / _delta(3, 2, vertices_omegas)) / 3;
}

static double _Ider_4(void)
{
  return 0.0;
}
