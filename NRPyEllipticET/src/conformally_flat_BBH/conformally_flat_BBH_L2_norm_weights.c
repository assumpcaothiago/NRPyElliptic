#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Faster replacement for repeated calls to conformally_flat_BBH_L2_norm_of_gf().
 *
 * That routine recomputes the radius and volume element (exp, sin, cos, sqrt) at every grid point on
 * every call; with info_output_freq = 0 it is called every relaxation step, where it costs a large
 * fraction of the step. The volume element and the r < integration_radius mask depend only on
 * (xx0, xx1), so we compute them once per grid, with exactly the same expressions, and reuse them.
 * The result is bitwise identical to conformally_flat_BBH_L2_norm_of_gf() (same terms, same order).
 */

// Returns malloc'ed weights w[i0 + Nxx0*i1] = dV (interior points with r < integration_radius), else 0.
REAL *conformally_flat_BBH_L2_norm_weights_malloc(const paramstruct *restrict params, const REAL integration_radius,
                                                  REAL *restrict xx[3]) {
#include "./conformally_flat_BBH_set_Cparameters.h"
  REAL *restrict w = (REAL *restrict)malloc(sizeof(REAL) * Nxx0 * Nxx1);
  if (w == NULL) return NULL;
  for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
    const REAL xx1 = xx[1][i1];
    for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
      const REAL xx0 = xx[0][i0];
      // Same expressions as conformally_flat_BBH_L2_norm_of_gf().
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
      w[(i0-NGHOSTS) + Nxx0*(i1-NGHOSTS)] = (r < integration_radius) ? sqrtdetgamma*dxx0*dxx1*dxx2 : 0.0;
    }
  }
  return w;
}

// log10 of the volume-weighted L2 norm of gridfunction gf_index, using precomputed weights.
REAL conformally_flat_BBH_L2_norm_of_gf_weighted(const paramstruct *restrict params, const int gf_index, const REAL *restrict w,
                                                 const REAL *restrict in_gf) {
#include "./conformally_flat_BBH_set_Cparameters.h"
  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;
  #pragma omp parallel for reduction(+:squared_sum,volume_sum)
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL dV = w[(i0-NGHOSTS) + Nxx0*(i1-NGHOSTS)];
        if (dV > 0.0) {
          const REAL gf_of_x = in_gf[IDX4S(gf_index,i0,i1,i2)];
          squared_sum += gf_of_x*gf_of_x*dV;
          volume_sum  += dV;
        }
      }
    }
  }
  return log10(1e-17 + sqrt(squared_sum/volume_sum));  // same formula as conformally_flat_BBH_L2_norm_of_gf()
}
