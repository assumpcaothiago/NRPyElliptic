#include "././NRPy_basic_defines.h"
/*
 * Reference Metric Precomputation infrastructure: Define rfmstruct
 */
void rfm_precompute_rfmstruct_define(const paramstruct *restrict params, REAL *restrict xx[3], rfm_struct *restrict rfmstruct) {
#include "./set_Cparameters.h"

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0[i0] = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))/(exp(1.0/SINHWAA) - exp(-1/SINHWAA));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__D0[i0] = AMAX*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)/(exp(1.0/SINHWAA) - exp(-1/SINHWAA));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__DD00[i0] = AMAX*(exp(xx0/SINHWAA)/pow(SINHWAA, 2) - exp(-xx0/SINHWAA)/pow(SINHWAA, 2))/(exp(1.0/SINHWAA) - exp(-1/SINHWAA));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__DDD000[i0] = AMAX*(exp(xx0/SINHWAA)/pow(SINHWAA, 3) + exp(-xx0/SINHWAA)/pow(SINHWAA, 3))/(exp(1.0/SINHWAA) - exp(-1/SINHWAA));
  }

  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct->f1_of_xx1[i1] = sin(xx1);
  }

  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct->f1_of_xx1__D1[i1] = cos(xx1);
  }

  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct->f1_of_xx1__DD11[i1] = -sin(xx1);
  }


  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct->f2_of_xx0_xx1[i0 + Nxx_plus_2NGHOSTS0*i1] = sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2));
  }


  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct->f2_of_xx0_xx1__D0[i0 + Nxx_plus_2NGHOSTS0*i1] = (1.0/2.0)*pow(AMAX, 2)*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2));
  }


  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct->f2_of_xx0_xx1__D1[i0 + Nxx_plus_2NGHOSTS0*i1] = pow(bScale, 2)*sin(xx1)*cos(xx1)/sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2));
  }


  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct->f2_of_xx0_xx1__DD00[i0 + Nxx_plus_2NGHOSTS0*i1] = -1.0/4.0*pow(AMAX, 4)*pow(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/(pow(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2), 3.0/2.0)*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 4)) + (1.0/2.0)*pow(AMAX, 2)*(2*exp(xx0/SINHWAA)/pow(SINHWAA, 2) - 2*exp(-xx0/SINHWAA)/pow(SINHWAA, 2))*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2)) + (1.0/2.0)*pow(AMAX, 2)*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2));
  }


  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct->f2_of_xx0_xx1__DD11[i0 + Nxx_plus_2NGHOSTS0*i1] = -pow(bScale, 4)*pow(sin(xx1), 2)*pow(cos(xx1), 2)/pow(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2), 3.0/2.0) - pow(bScale, 2)*pow(sin(xx1), 2)/sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2)) + pow(bScale, 2)*pow(cos(xx1), 2)/sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2)*pow(sin(xx1), 2));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f3_of_xx0[i0] = sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f3_of_xx0__D0[i0] = (1.0/2.0)*pow(AMAX, 2)*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2));
  }

  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f3_of_xx0__DD00[i0] = -1.0/4.0*pow(AMAX, 4)*pow(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/(pow(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2), 3.0/2.0)*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 4)) + (1.0/2.0)*pow(AMAX, 2)*(2*exp(xx0/SINHWAA)/pow(SINHWAA, 2) - 2*exp(-xx0/SINHWAA)/pow(SINHWAA, 2))*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2)) + (1.0/2.0)*pow(AMAX, 2)*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)/(sqrt(pow(AMAX, 2)*pow(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA), 2)/pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2) + pow(bScale, 2))*pow(exp(1.0/SINHWAA) - exp(-1/SINHWAA), 2));
  }

}
