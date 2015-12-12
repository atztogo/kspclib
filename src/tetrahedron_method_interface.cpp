#include <iostream>
#include <vector>
extern "C" {
#include "kgrid.h"
#include "tetrahedron_method.h"
}

static int grid_address_to_index(int [3], int [3]);

void calc_dos_lin_tetra(double *energies, int *grid_address, int *spg_grid_mapping_table, int *mesh, double *rec_basis, double *energy_samples, int num_energy_samples,int num_bands, int num_kpoints_ibz, double volume, int bloechl, double *dos, double *int_dos);

void calc_dos_lin_tetra(double *energies, int *grid_address, int *spg_grid_mapping_table, int *mesh, double *rec_basis, double *energy_samples, int num_energy_samples,int num_bands, int num_kpoints_ibz, double volume, int bloechl, double *dos, double *int_dos) {

  int num_kpoints=mesh[0]*mesh[1]*mesh[2];
  int weights[num_kpoints];
  int ir_gp[num_kpoints_ibz];
  std::vector<int> ir_weights(num_kpoints_ibz);
  int gp_ir_index[num_kpoints];
  double omegas[24][4];
  int g_addr[3];
  int relative_grid_address[24][4][3];

  // spglib does not use pointers or vectors, so convert 2D vectors to std. c arrays
  // (1d vector is compatible with a 1d array, so leave those intact)                        
  std::vector<std::vector<double> > energies_temp(num_bands, std::vector<double> (num_kpoints_ibz));
  for(int band=0;band<num_bands;band++) {
    for (int kpoint=0;kpoint<num_kpoints_ibz;kpoint++) {
      energies_temp[band][kpoint] = *((double *)energies + band*num_kpoints_ibz + kpoint);
    }
  }
  double rec_basis_temp[3][3];
  for(int i=0;i<3;++i) {
    for (int j=0;j<3;j++) {
      rec_basis_temp[i][j] = *((double *)rec_basis + i*3 + j);
    }
  }
  double grid_address_temp[num_kpoints][3];
  for(int kpoint=0;kpoint<num_kpoints;kpoint++) {
    for (int dir=0;dir<3;dir++) {
      grid_address_temp[kpoint][dir] = *((int *)grid_address + kpoint*3 + dir);
    }
  }

  // initialize weights for BZ
  for (int kpoint = 0; kpoint < num_kpoints; kpoint++) {
    weights[kpoint] = 0.0;
  }

  // set weights for BZ
  for (int kpoint = 0; kpoint < num_kpoints; kpoint++) {
    weights[spg_grid_mapping_table[kpoint]]++;
  }

  // set weights in IBZ and the kmesh mapping indexes
  for (int kpoint = 0, index=0; kpoint < num_kpoints; kpoint++) {
    if (weights[kpoint] != 0) {
      ir_gp[index] = kpoint;
      ir_weights[index] = weights[kpoint];
      gp_ir_index[kpoint] = index;
      index++;
    } 
    else {
      gp_ir_index[kpoint] = gp_ir_index[spg_grid_mapping_table[kpoint]];
    }
  }
  thm_get_relative_grid_address(relative_grid_address, rec_basis_temp);

  // now find the integration weights
  for (int energy = 0; energy < num_energy_samples; energy++) {
    for (int band=0; band<num_bands; band++) {
      dos[band*num_energy_samples+energy]=0.0;
      int_dos[band*num_energy_samples+energy]=0.0;
      for (int kpoint = 0; kpoint < num_kpoints_ibz;  kpoint++) {
        for (int tetra = 0; tetra < 24; tetra++) {
          for (int corner = 0; corner < 4; corner++) {
            for (int dir = 0; dir < 3; dir++) {
              g_addr[dir] = grid_address_temp[ir_gp[kpoint]][dir] +
                relative_grid_address[tetra][corner][dir];
            }
            int gp = grid_address_to_index(g_addr, mesh);
	    omegas[tetra][corner] = energies_temp[band][gp_ir_index[gp]];
          }
        }
        dos[band*num_energy_samples+energy] += ir_weights[kpoint] * thm_get_integration_weight(energy_samples[energy], omegas, 'I', bloechl);
        int_dos[band*num_energy_samples+energy] += ir_weights[kpoint] * thm_get_integration_weight(energy_samples[energy], omegas, 'J', bloechl);
      }
    } 
  }
}

static int grid_address_to_index(int g[3], int mesh[3]) {
  int i;
  int gm[3];

  for (i = 0; i < 3; i++) {
    gm[i] = g[i] % mesh[i];
    if (gm[i] < 0) {
      gm[i] += mesh[i];
    }
  }
  return (gm[0] + gm[1] * mesh[0] + gm[2] * mesh[0] * mesh[1]);
}
