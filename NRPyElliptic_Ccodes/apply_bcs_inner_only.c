#include "././NRPy_basic_defines.h"
/*
 * 
 * Apply BCs to inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 */
void apply_bcs_inner_only(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs) {
#include "./set_Cparameters.h"


  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

#pragma omp parallel for  // spawn threads and distribute across them
  for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;

      gfs[IDX4ptS(which_gf, dstpt)] = bcstruct->inner_bc_array[pt].parity[evol_gf_parity[which_gf]] * gfs[IDX4ptS(which_gf, srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++)
}
