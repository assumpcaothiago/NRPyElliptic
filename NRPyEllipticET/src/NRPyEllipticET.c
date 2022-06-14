// NRPyEllipticET header files
#include "NRPyEllipticET.h"

CCTK_REAL *NRPyEllipticET_uu = NULL;
CCTK_REAL *NRPyEllipticET_xx[3] = {NULL,NULL,NULL};

void NRPyEllipticET(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Set the size of the local grid
  const CCTK_INT Nxx_plus_2NGHOSTS0    = N0 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS1    = N1 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS2    = N2 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  if( NRPyEllipticET_uu == NULL ) {
    // The following instructions only happen the
    // first time this function is called.

    // Allocate memory for NRPyEllipticET_xx and NRPyEllipticET_uu
    NRPyEllipticET_xx[0] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS0);
    NRPyEllipticET_xx[1] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS1);
    NRPyEllipticET_xx[2] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS2);
    NRPyEllipticET_uu    = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS_tot);

    // Check everything is okay
    if( !NRPyEllipticET_uu || !NRPyEllipticET_xx[0] || !NRPyEllipticET_xx[1] || !NRPyEllipticET_xx[2] ) {
      CCTK_ERROR("Could not allocate memory for NRPyEllipticET_uu or NRPyEllipticET_xx!");
    }

    // Now solve the constraints of General Relativity
    NRPyEllipticET_PunctureInitialData_Hyperbolic_Relaxation();
  }

  // Check memory has been allocated for NRPyEllipticET_uu. Error out otherwise.
  if( !NRPyEllipticET_uu || !NRPyEllipticET_xx[0] || !NRPyEllipticET_xx[1] || !NRPyEllipticET_xx[2] ) {
    CCTK_ERROR("Reached NRPyEllipticET_PunctureInitialData_Initialize_ADMBase() but memory for NRPyEllipticET_uu or NRPyEllipticET_xx has not been allocated.");
  }
  else {
    // Now interpolate the solution to the local ETK grid
    NRPyEllipticET_PunctureInitialData_Initialize_ADMBase(cctkGH,cctk_lsh,x,y,z,uuGF,
                                                          alp,betax,betay,betaz,
                                                          gxx,gxy,gxz,gyy,gyz,gzz,
                                                          kxx,kxy,kxz,kyy,kyz,kzz);
  }
}
