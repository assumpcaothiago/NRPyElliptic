#include "./conformally_flat_BBH_NRPy_basic_defines.h"
/*
 * Find the CFL-constrained timestep
 */
REAL conformally_flat_BBH_find_timestep(const paramstruct *restrict params, REAL *restrict xx[3], const REAL CFL_FACTOR) {
#include "./conformally_flat_BBH_set_Cparameters.h"
    REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL xx0 = xx[0][i0];
        REAL ds_dirn0, ds_dirn1, ds_dirn2;
        /*
         *  Original SymPy expressions:
         *  "[ds_dirn0 = AMAX*dxx0*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))),
         *    ds_dirn1 = dxx1*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2),
         *    ds_dirn2 = AMAX*dxx2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))]"
         */
        {
          const double tmp_1 = sin(xx1);
          const double tmp_2 = (1.0/(SINHWAA));
          const double tmp_3 = exp(tmp_2) - exp(-tmp_2);
          const double tmp_5 = exp(tmp_2*xx0);
          const double tmp_6 = exp(-tmp_2*xx0);
          const double tmp_7 = tmp_5 - tmp_6;
          const double tmp_8 = ((AMAX)*(AMAX))*((tmp_7)*(tmp_7))/((tmp_3)*(tmp_3));
          const double tmp_9 = sqrt(((bScale)*(bScale))*((tmp_1)*(tmp_1)) + tmp_8);
          const double tmp_10 = AMAX/tmp_3;
          ds_dirn0 = dxx0*tmp_10*tmp_9*(tmp_2*tmp_5 + tmp_2*tmp_6)/sqrt(((bScale)*(bScale)) + tmp_8);
          ds_dirn1 = dxx1*tmp_9;
          ds_dirn2 = dxx2*tmp_1*tmp_10*tmp_7;
        }

        #ifndef MIN
        #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
        #endif
        // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2):
        dsmin = MIN(dsmin, MIN(ds_dirn0, MIN(ds_dirn1, ds_dirn2)));

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    return dsmin/1.0;
}
