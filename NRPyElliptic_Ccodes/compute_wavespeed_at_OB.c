#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Compute wavespeed at outer boundary
 */
REAL compute_wavespeed_at_OB(const paramstruct *restrict params, REAL *restrict auxevol_gfs) {
#include "./set_Cparameters.h"


  REAL wavespeed_at_OB = auxevol_gfs[IDX4S(WAVESPEEDGF, Nxx_plus_2NGHOSTS0-NGHOSTS-1, NGHOSTS, Nxx_plus_2NGHOSTS2/2)];
  return wavespeed_at_OB;
}
