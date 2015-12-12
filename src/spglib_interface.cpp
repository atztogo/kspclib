#include <iostream>
extern "C" {
#include <spglib.h>
}


void get_reciprocal_mesh_interface(int *mesh, double *lattice, double *positions, int *anumbers, unsigned int num_atoms, int *is_shifted, int *mesh_points, int *mapping, int is_time_reversal, double symprec);

void get_reciprocal_mesh_interface(int *mesh, double *lattice, double *positions, int *anumbers, unsigned int num_atoms, int *is_shifted, int *mesh_points, int *mapping, int is_time_reversal, double symprec) {

  const unsigned int num_k_points=mesh[0]*mesh[1]*mesh[2];
  unsigned int i,j;

  // spglib does not use pointers or vectors, so convert 2D vectors to std. c arrays
  // (1d vector is compatible with a 1d array, so leave those intact)
  int temp_mesh_points[num_k_points][3];
  for(i=0;i<num_k_points;++i) {
    for (j=0;j<3;j++) {
      temp_mesh_points[i][j] = *((int *)mesh_points + i*3 + j);
    }
  }
  double temp_lattice[3][3];
  for(i=0;i<3;++i) {
    for (j=0;j<3;j++) {
      temp_lattice[i][j]=*((double *)lattice + i*3 + j);
    }
  }
  double temp_positions[num_atoms][3];
  for(i=0;i<num_atoms;++i) {
    for (j=0;j<3;j++) {
      temp_positions[i][j]=*((double *)positions + i*3 + j);
    }
  }
  // then get the kpoint mesh and associated ibz through the mapping table
  spg_get_ir_reciprocal_mesh(temp_mesh_points,mapping,mesh,is_shifted,is_time_reversal,temp_lattice,temp_positions,anumbers,num_atoms,symprec);
  for(i=0;i<num_k_points;++i) {
    for (j=0;j<3;j++) {
      *((int *)mesh_points + i*3 + j)=temp_mesh_points[i][j];
    }
  }
}
