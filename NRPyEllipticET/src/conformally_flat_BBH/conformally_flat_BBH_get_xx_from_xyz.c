#include "./conformally_flat_BBH_NRPy_basic_defines.h"

/*
 * This function uses a Newton-Raphson solver to find (xx0,xx1,xx2) from (x,y,z)
 */
void conformally_flat_BBH_get_xx_from_xyz(const double AMAX, const double SINHWAA, const double bScale,
                                          const double xCart[3], double xx[3]) {

      // Set local Cartesian variables
      const double Cartx = xCart[0];
      const double Carty = xCart[1];
      const double Cartz = xCart[2];
  const double tmp_1 = ((bScale)*(bScale));
  const double tmp_2 = ((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz));
  const double tmp_3 = sqrt(-4*((Cartz)*(Cartz))*tmp_1 + ((bScale)*(bScale)*(bScale)*(bScale)) + 2*tmp_1*tmp_2 + ((tmp_2)*(tmp_2)));
  const double tmp_4 = (1.0/(tmp_1));
  xx[0] = SINHWAA*asinh(M_SQRT1_2*sqrt(-tmp_1 + tmp_2 + tmp_3)*sinh((1.0/(SINHWAA)))/AMAX);
  xx[1] = acos(0.99999999*M_SQRT1_2*sqrt(tmp_2*tmp_4 - tmp_3*tmp_4 + 1)*(((Cartz) > 0) - ((Cartz) < 0)));
  xx[2] = atan2(Carty, Cartx);
}
