// NRPyEllipticET header files
#include "NRPyEllipticET.h"

CCTK_REAL *NRPyEllipticET_uu = NULL;
CCTK_REAL *NRPyEllipticET_xx[3] = {NULL,NULL,NULL};

void NRPyEllipticET(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We only find the solution to the elliptic
  // problem on the first time this function is
  // called (typically on the coarsest Carpet
  // refinement level). The following if statement
  // takes care of that. Note, however, that the
  // solution will be obtained once per MPI process.
  if( NRPyEllipticET_uu == NULL ) {
    // Now solve the constraints of General Relativity
    if( CCTK_Equals(initial_data_type, "ConformallyFlatBBH") ) {
      NRPyEllipticET_conformally_flat_BBH_Hyperbolic_Relaxation();
    }
  }

  // Check memory has been allocated for NRPyEllipticET_uu. Error out otherwise.
  if( !NRPyEllipticET_uu || !NRPyEllipticET_xx[0] || !NRPyEllipticET_xx[1] || !NRPyEllipticET_xx[2] ) {
    CCTK_ERROR("Reached NRPyEllipticET_PunctureInitialData_Initialize_ADMBase() but memory for NRPyEllipticET_uu or NRPyEllipticET_xx has not been allocated.");
  }
  else {
    // Now interpolate the solution to the local ETK grid
    if( CCTK_Equals(initial_data_type, "ConformallyFlatBBH") ) {
      NRPyEllipticET_conformally_flat_BBH_Initialize_ADMBase(cctkGH,cctk_lsh,x,y,z,uuGF,
                                                             alp,betax,betay,betaz,
                                                             gxx,gxy,gxz,gyy,gyz,gzz,
                                                             kxx,kxy,kxz,kyy,kyz,kzz);
    }
  }
}
