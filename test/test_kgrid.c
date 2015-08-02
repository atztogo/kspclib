/* Function: kgd_get_all_grid_addresses */
/*           Show indeces of grid points and their addresses */
/* */
/* Compile1: Run faster from left most index */
/*     % gcc test_kgrid.c [source_dir]/kgrid.c -I[source_dir] */
/* Compile2: Run faster from left most index */
/*     % gcc -DGRID_ORDER_XYZ test_kgrid.c [source_dir]/kgrid.c -I[source_dir] */

#include "kgrid.h"
#include <stdio.h>

int main(void)
{
  int i;
  int mesh[3] = {5, 5, 5};
  int is_shift[3] = {1, 1, 1};
  int address_double[3];
  int num_gp = mesh[0] * mesh[1] * mesh[2];
  int grid_address[num_gp][3];

  kgd_get_all_grid_addresses(grid_address, mesh);

  printf("index :    grid     : double-grid (shift)\n");
  for (i = 0; i < num_gp; i ++) {
    kgd_get_grid_address_double_mesh(address_double,
				     grid_address[i],
				     mesh,
				     is_shift);
    printf("%5d : %3d %3d %3d : %3d %3d %3d (%d %d %d)\n", i,
	   grid_address[i][0],
	   grid_address[i][1],
	   grid_address[i][2],
	   address_double[0],
	   address_double[1],
	   address_double[2],
	   (is_shift[0] != 0),
	   (is_shift[1] != 0),
	   (is_shift[2] != 0));
  }


  return 0;
}
