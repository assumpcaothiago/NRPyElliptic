#ifndef NRPYELLIPTICET_H_
#define NRPYELLIPTICET_H_

// NRPy Interface
#include "NRPy_basic_defines.h"
#include "NRPy_function_prototypes.h"

// Standard C includes
#include <time.h>

// This pointer will contain the solution for u
//__attribute__ ((unused)) static CCTK_REAL* NRPyEllipticET_uu = NULL;
extern CCTK_REAL *NRPyEllipticET_uu;

// This pointer contains the numerical grid used by NRPyEllipticET
//__attribute__ ((unused)) static CCTK_REAL* NRPyEllipticET_xx[3] = {NULL,NULL,NULL};
extern CCTK_REAL *NRPyEllipticET_xx[3];

#endif // NRPYELLIPTICET_H_
