#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Curvilinear boundary condition driver routine: Apply BCs to all six
 *   boundary faces of the 3D numerical domain, filling in the
 *   innermost ghost zone layer first, and moving outward.
 */
void apply_bcs_curvilinear_inner_only(const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict xx[3],
                           REAL *restrict gfs) {
#include "./set_Cparameters.h"

#pragma omp parallel for
  for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      // Apply INNER (parity) boundary conditions:
      for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
        const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
        const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
        const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
        const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
        const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
        const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
        gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
          bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src,i1src,i2src)];
      } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
    } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
  } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
}
