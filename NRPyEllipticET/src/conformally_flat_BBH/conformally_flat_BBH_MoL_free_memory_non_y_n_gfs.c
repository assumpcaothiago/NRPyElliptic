#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Method of Lines (MoL) for "RK4" method: Free memory for "non_y_n_gfs" gridfunctions
 *    - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
 *    - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
 */
void conformally_flat_BBH_MoL_free_memory_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs) {
#include "./conformally_flat_BBH_set_Cparameters.h"

      free(gridfuncs->y_nplus1_running_total_gfs);
      free(gridfuncs->k_odd_gfs);
      free(gridfuncs->k_even_gfs);
      free(gridfuncs->auxevol_gfs);
}
