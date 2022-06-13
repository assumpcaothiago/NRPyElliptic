#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Update evolgf_speed at outer boundary
 */
void update_evolgf_speed(const paramstruct *restrict params, REAL *restrict auxevol_gfs, REAL *evolgf_speed) {
#include "./set_Cparameters.h"


  REAL wavespeed_at_OB = auxevol_gfs[IDX4S(WAVESPEEDGF, Nxx_plus_2NGHOSTS0 - NGHOSTS - 1, 
                                     NGHOSTS, Nxx_plus_2NGHOSTS2/2)];
  for (int ii = 0; ii < NUM_EVOL_GFS; ii++){
    evolgf_speed[ii] = wavespeed_at_OB;
  } // END for
}
