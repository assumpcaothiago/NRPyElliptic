#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Free memory allocated within bcstruct
 */
void freemem_bcstruct(const paramstruct *restrict params, const bc_struct *restrict bcstruct) {
#include "./set_Cparameters.h"

  for(int i=0;i<NGHOSTS;i++) { free(bcstruct->outer[i]);  free(bcstruct->inner[i]); }
  free(bcstruct->outer);  free(bcstruct->inner);
  free(bcstruct->num_ob_gz_pts); free(bcstruct->num_ib_gz_pts);
}
