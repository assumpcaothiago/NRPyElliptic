#include "././NRPy_basic_defines.h"
/*
 * Reference Metric Precomputation infrastructure: Free rfmstruct memory
 */
void rfm_precompute_rfmstruct_freemem(const paramstruct *restrict params, rfm_struct *restrict rfmstruct) {
#include "./set_Cparameters.h"

  free(rfmstruct->f0_of_xx0);
  free(rfmstruct->f0_of_xx0__D0);
  free(rfmstruct->f0_of_xx0__DD00);
  free(rfmstruct->f0_of_xx0__DDD000);
  free(rfmstruct->f1_of_xx1);
  free(rfmstruct->f1_of_xx1__D1);
  free(rfmstruct->f1_of_xx1__DD11);
  free(rfmstruct->f2_of_xx0_xx1);
  free(rfmstruct->f2_of_xx0_xx1__D0);
  free(rfmstruct->f2_of_xx0_xx1__D1);
  free(rfmstruct->f2_of_xx0_xx1__DD00);
  free(rfmstruct->f2_of_xx0_xx1__DD11);
  free(rfmstruct->f3_of_xx0);
  free(rfmstruct->f3_of_xx0__D0);
  free(rfmstruct->f3_of_xx0__DD00);
}
