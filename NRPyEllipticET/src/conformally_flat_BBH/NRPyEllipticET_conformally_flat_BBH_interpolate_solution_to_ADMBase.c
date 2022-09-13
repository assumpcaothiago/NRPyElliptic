#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "util_Table.h"

void NRPyEllipticET_conformally_flat_BBH_interpolate_solution_to_ADMBase( const CCTK_INT num_pts,
                                                                          const CCTK_INT *input_array_dims,
                                                                          const CCTK_REAL *origin,
                                                                          const CCTK_REAL *deltas,
                                                                          const void **interp_coords,
                                                                          const CCTK_REAL *input_gf,
                                                                          CCTK_REAL *output_gf ) {

  DECLARE_CCTK_PARAMETERS;

  // Get interpolation operator handle
  CCTK_INT operator_handle = CCTK_InterpHandle(interpolator_name);
  if( operator_handle < 0 ) CCTK_WARN(CCTK_WARN_ABORT,"Error while getting interpolation operator handle.");

  // Set interpolation order
  char interp_order_string[10];
  sprintf(interp_order_string,"order=%d",interpolation_order);
  CCTK_INT param_table_handle = Util_TableCreateFromString(interp_order_string);
  if( param_table_handle < 0 ) CCTK_WARN(CCTK_WARN_ABORT,"Error while creating parameter table.");

  // Set input array variable types
  CCTK_INT input_array_type_codes[1] = { CCTK_VARIABLE_REAL };

  // Set output array variable types
  CCTK_INT output_array_type_codes[1] = { CCTK_VARIABLE_REAL };

  // Set pointer to input array
  const void* input_arrays[1] = {
    (const void *)input_gf
  };

  // Set output gridfunctions
  void* output_arrays[1] = {
    (void *)output_gf
  };

  CCTK_INT ierr = CCTK_InterpLocalUniform(3,
                                          operator_handle,param_table_handle,
                                          origin,deltas,
                                          num_pts,
                                          CCTK_VARIABLE_REAL,
                                          interp_coords,
                                          1,
                                          input_array_dims,
                                          input_array_type_codes,
                                          input_arrays,
                                          1,
                                          output_array_type_codes,
                                          output_arrays);

  if( ierr < 0 ) {
    Util_TableDestroy(param_table_handle);
    CCTK_VWARN(CCTK_WARN_ABORT,"Interpolation screwed up.");
  }
  else {
    ierr = Util_TableDestroy(param_table_handle);
    if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT,"Could not destroy table");
  }
}
