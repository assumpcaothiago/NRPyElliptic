#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Free memory allocated within bcstruct
 */
void conformally_flat_BBH_freemem_bcstruct(const paramstruct *restrict params, const bc_struct *restrict bcstruct) {
#include "./conformally_flat_BBH_set_Cparameters.h"

  for(int i=0;i<NGHOSTS;i++) { free(bcstruct->outer[i]);  free(bcstruct->inner[i]); }
  free(bcstruct->outer);  free(bcstruct->inner);
  free(bcstruct->num_ob_gz_pts); free(bcstruct->num_ib_gz_pts);
}
