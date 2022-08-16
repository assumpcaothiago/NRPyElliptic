#include "./conformally_flat_BBH_NRPy_basic_defines.h"

CCTK_REAL conformally_flat_BBH_compute_TRC(const CCTK_REAL A,
                                           const CCTK_REAL b,
                                           const CCTK_REAL C_0,
                                           const CCTK_REAL dt,
                                           const CCTK_REAL dxx1,
                                           const CCTK_REAL dxx2) {
  // Eq. (B6) of https://arxiv.org/pdf/2111.02424.pdf
  // Note that we use (xx0,xx1,xx2) while the paper
  // uses (x1,x2,x3).
  const CCTK_REAL dx1o2    = dxx1/2.0;
  const CCTK_REAL dx2o2    = dxx2/2.0;
  const CCTK_REAL cosdx1o2 = cos(dx1o2);
  const CCTK_REAL tmp1     = A / sqrt( A*A + b*b*cosdx1o2*cosdx1o2 );
  return dt / (C_0 * dxx1) * atanh(tmp1) * cos(dx1o2) * cos(dx2o2);
}

void conformally_flat_BBH_set_optimum_eta_damping(const CCTK_REAL dt,
                                                  const CCTK_REAL CFL_FACTOR,
                                                  paramstruct *restrict params) {
#include "conformally_flat_BBH_set_Cparameters.h"

  const CCTK_REAL TRC = conformally_flat_BBH_compute_TRC(AMAX, bScale, CFL_FACTOR, dt, dxx1, dxx2);
  // Eq. (36) of https://arxiv.org/pdf/2111.02424.pdf
  params->eta_damping = (6.9*pow(TRC,-0.79) - 0.7);
}
