#include "././NRPy_basic_defines.h"
/*
 * Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2),   accounting for the origin of this grid being possibly offcenter.
 */
void xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include "./set_Cparameters.h"


    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
      /*
   *  Original SymPy expressions:
   *  "[xCart[0] = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*cos(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)) + Cart_originx,
   *    xCart[1] = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*sin(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)) + Cart_originy,
   *    xCart[2] = Cart_originz + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)]"
   */
  {
    const double tmp_0 = (1.0/(SINHWAA));
    const double tmp_1 = exp(tmp_0) - exp(-tmp_0);
    const double tmp_3 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
    const double tmp_4 = AMAX*tmp_3*sin(xx1)/tmp_1;
    xCart[0] = Cart_originx + tmp_4*cos(xx2);
    xCart[1] = Cart_originy + tmp_4*sin(xx2);
    xCart[2] = Cart_originz + sqrt(((AMAX)*(AMAX))*((tmp_3)*(tmp_3))/((tmp_1)*(tmp_1)) + ((bScale)*(bScale)))*cos(xx1);
  }
}
