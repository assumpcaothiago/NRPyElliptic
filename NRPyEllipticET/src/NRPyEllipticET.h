#ifndef NRPYELLIPTICET_H_
#define NRPYELLIPTICET_H_

// Standard C includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Einstein Toolkit interface
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// This pointer will contain the solution for u
extern CCTK_REAL *NRPyEllipticET_uu;

// This pointer contains the numerical grid used by NRPyEllipticET
extern CCTK_REAL *NRPyEllipticET_xx[3];

// Function prototypes
// .-------------------------------------------------.
// | Conformally flat binary black hole initial data |
// .-------------------------------------------------.
void NRPyEllipticET_conformally_flat_BBH_Hyperbolic_Relaxation();
void NRPyEllipticET_conformally_flat_BBH_Initialize_ADMBase( const cGH *restrict cctkGH,
                                                             const int *restrict cctk_lsh,
                                                             const CCTK_REAL *restrict x,
                                                             const CCTK_REAL *restrict y,
                                                             const CCTK_REAL *restrict z,
                                                             CCTK_REAL *restrict uuGF,
                                                             CCTK_REAL *restrict alp,
                                                             CCTK_REAL *restrict betax,
                                                             CCTK_REAL *restrict betay,
                                                             CCTK_REAL *restrict betaz,
                                                             CCTK_REAL *restrict gxx,
                                                             CCTK_REAL *restrict gxy,
                                                             CCTK_REAL *restrict gxz,
                                                             CCTK_REAL *restrict gyy,
                                                             CCTK_REAL *restrict gyz,
                                                             CCTK_REAL *restrict gzz,
                                                             CCTK_REAL *restrict kxx,
                                                             CCTK_REAL *restrict kxy,
                                                             CCTK_REAL *restrict kxz,
                                                             CCTK_REAL *restrict kyy,
                                                             CCTK_REAL *restrict kyz,
                                                             CCTK_REAL *restrict kzz );

#endif // NRPYELLIPTICET_H_
