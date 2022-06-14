#include "NRPyEllipticET.h"

void NRPyEllipticET_cleanup(CCTK_ARGUMENTS) {
  // Free memory allocated for NRPyEllipticET_uu and NRPyEllipticET_xx
  if( NRPyEllipticET_uu != NULL ) {
    free(NRPyEllipticET_uu); NRPyEllipticET_uu = NULL;
    for(int i=0;i<3;i++) {
      free(NRPyEllipticET_xx[i]); NRPyEllipticET_xx[i] = NULL;
    }
  }
  // Check everything is okay
  if( (NRPyEllipticET_uu    == NULL) &&
      (NRPyEllipticET_xx[0] == NULL) &&
      (NRPyEllipticET_xx[1] == NULL) &&
      (NRPyEllipticET_xx[2] == NULL) ) {
    CCTK_INFO("Successfully freed memory for NRPyEllipticET_uu and NRPyEllipticET_xx!");
  }
  else {
    CCTK_ERROR("Problem freeing memory for NRPyEllipticET_uu or NRPyEllipticET_xx!");
  }
  // All done!
  CCTK_INFO("All done!");
}
