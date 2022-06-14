#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * driver_bcstruct(): Set up bcstruct.
 * WARNING: Needs Nxx_plus_2NGHOSTS_tot * [sizeof(gz_map) + sizeof(parity_condition) ] bytes of temporary memory!
 *    (This memory will be deallocated at end of function, in freemem_bcstruct().)
 *    Note that gz_map consists of 3*sizeof(short = 2 bytes), and parity_condition consists of 10*sizeof(int8_t = 1 byte)
 *    Thus the total cost is about 16 bytes at all gridpoints, or 2 double-precision gridfunctions.
 *    To avoid memory overload, be sure to set this up prior to allocation of other gridfunctions.
 * STEPS:
 * 1. Allocate memory for bc_gz_map, which maps inner gzs to the
 *    appropriate interior gridpoint. In the case of outer gzs, the
 *    gz point maps to itself.
 * 2. Allocate storage for parity_condition, which
 */
void driver_bcstruct(const paramstruct *restrict params, bc_struct *restrict bcstruct, REAL *restrict xx[3]) {
#include "./set_Cparameters.h"


  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  // Step 1: Allocate memory storage for bc_gz_map, which
  //         in the case a boundary point is a *parity*
  //         boundary, is set to the interior, non-
  //         boundary point corresponding to the same
  //         Cartesian gridpoint. Otherwise bc_gz_map
  //         is set to (i0,i1,i2) = (-1,-1,-1).
  gz_map *restrict bc_gz_map = (gz_map *restrict)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);

  // Step 2: Allocate memory storage for bc_parity_conditions,
  //         which store parity conditions for all 10
  //         gridfunction types at all grid points.
  parity_condition *restrict bc_parity_conditions = (parity_condition *restrict)malloc(sizeof(parity_condition)*Nxx_plus_2NGHOSTS_tot);

  // Step 3: Set bc_gz_map and bc_parity_conditions at *all*
  //         points; on the boundary and otherwise.
  set_up__bc_gz_map_and_parity_condns(params, xx, bc_gz_map,
                                      bc_parity_conditions);

  // Step 4: Declare and allocate memory for bcstruct,
  //         which will store all information needed for
  //         applying the boundary conditions.
  bcstruct->outer = (outer_bc **)malloc(sizeof(outer_bc *)*NGHOSTS);
  bcstruct->inner = (inner_bc **)malloc(sizeof(inner_bc *)*NGHOSTS);
  bcstruct->num_ob_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);
  bcstruct->num_ib_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);

  // Step 5: Store all information needed to quickly and
  //         efficiently apply boundary conditions. This
  //         function transfers all information from
  //         bc_gz_map (defined at *all gridpoints*) into
  //         bcstruct (defined only at boundary points).
  //         Thus when this function has finished,
  //         bc_gz_map is no longer needed.
  set_bcstruct(params,bc_gz_map,
               bc_parity_conditions,
               bcstruct);
  // Do not apply outer boundary conditions on grids that do not possess the outer boundary!
  //if(params->has_outer_boundary == 0) {
  //  for(int i=0;i<NGHOSTS;i++) bcstruct->num_ob_gz_pts[i] = 0;
  //}

  // Step 6: As described in Step 4, bc_gz_map is no
  //         longer needed at this point, so we free its
  //         memory. Farewell, friend!
  free(bc_gz_map);
  free(bc_parity_conditions);
}
