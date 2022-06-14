#include "././NRPy_basic_defines.h"
/*
 * Reference Metric Precomputation infrastructure: Allocate memory for rfmstruct
 */
void rfm_precompute_rfmstruct_malloc(const paramstruct *restrict params, rfm_struct *restrict rfmstruct) {
#include "./set_Cparameters.h"

  rfmstruct->f0_of_xx0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__D0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DDD000 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f1_of_xx1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__D1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__DD11 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
  rfmstruct->f2_of_xx0_xx1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1);
  rfmstruct->f2_of_xx0_xx1__D0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1);
  rfmstruct->f2_of_xx0_xx1__D1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1);
  rfmstruct->f2_of_xx0_xx1__DD00 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1);
  rfmstruct->f2_of_xx0_xx1__DD11 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1);
  rfmstruct->f3_of_xx0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f3_of_xx0__D0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
  rfmstruct->f3_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
}
