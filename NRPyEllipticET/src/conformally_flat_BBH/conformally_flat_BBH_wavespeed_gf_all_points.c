#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Maximum CFL-constrained wavespeed at all points
 */
void conformally_flat_BBH_wavespeed_gf_all_points(const paramstruct *restrict params, const REAL CFL_FACTOR, const REAL dt, REAL *restrict xx[3], REAL *restrict auxevol_gfs) {
#include "./conformally_flat_BBH_set_Cparameters.h"

  #pragma omp parallel for
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

        #ifndef MIN
        #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
        #endif
        // Set dsmin = minimum proper distance;
        REAL dsmin = MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2));
        auxevol_gfs[IDX4S(WAVESPEEDGF, i0, i1, i2)] = CFL_FACTOR*dsmin/dt;
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
