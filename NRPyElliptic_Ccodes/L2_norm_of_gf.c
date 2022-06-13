#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Evaluate L2-norm of gridfunction
 */
REAL L2_norm_of_gf(const paramstruct *restrict params, const int gf_index, const REAL integration_radius, REAL *restrict xx[3], const REAL *restrict in_gf) {
#include "./set_Cparameters.h"

  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;
    
  #pragma omp parallel for reduction(+:squared_sum,volume_sum)
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL xx0 = xx[0][i0];
        const double tmp_0 = sin(xx1);
        const double tmp_3 = (1.0/(SINHWAA));
        const double tmp_4 = exp(tmp_3) - exp(-tmp_3);
        const double tmp_6 = exp(tmp_3*xx0);
        const double tmp_7 = exp(-tmp_3*xx0);
        const double tmp_8 = ((tmp_6 - tmp_7)*(tmp_6 - tmp_7));
        const double tmp_9 = ((AMAX)*(AMAX))*tmp_8/((tmp_4)*(tmp_4));
        const double tmp_11 = ((bScale)*(bScale)) + tmp_9;
        const REAL r = sqrt(((tmp_0)*(tmp_0))*tmp_9 + tmp_11*((cos(xx1))*(cos(xx1))));
        const REAL sqrtdetgamma = ((AMAX)*(AMAX))*sqrt(tmp_8*((((bScale)*(bScale))*((tmp_0)*(tmp_0)) + tmp_9)*(((bScale)*(bScale))*((tmp_0)*(tmp_0)) + tmp_9))*((tmp_3*tmp_6 + tmp_3*tmp_7)*(tmp_3*tmp_6 + tmp_3*tmp_7))/(tmp_11*((tmp_4)*(tmp_4)*(tmp_4)*(tmp_4))))*fabs(tmp_0);
        
        if(r < integration_radius) {
          const REAL gf_of_x = in_gf[IDX4S(gf_index,i0,i1,i2)];
          const REAL dV = sqrtdetgamma*dxx0*dxx1*dxx2;
          squared_sum += gf_of_x*gf_of_x*dV;
          volume_sum  += dV;
        } // END if(r < integration_radius)
        
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)

  // Compute and output the log of the L2 norm.
  return log10(1e-16 + sqrt(squared_sum/volume_sum));  // 1e-16 + ... avoids log10(0)
}
