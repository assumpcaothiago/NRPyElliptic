#include "././NRPy_basic_defines.h"
/*
 * Given Cartesian point (x,y,z), this function outputs the corresponding
 *   (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2_global_grid_center(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]) {
#include "./set_Cparameters.h"

    const REAL Cartx = xCart[0];
    const REAL Carty = xCart[1];
    const REAL Cartz = xCart[2];

  /*
   *  Original SymPy expressions:
   *  "[xx[0] = M_SQRT1_2*sqrt(Cartx**2 + Carty**2 + Cartz**2 - bScale**2 + sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)),
   *    xx[1] = acos(M_SQRT1_2*sqrt(1 + (Cartx**2 + Carty**2 + Cartz**2)/bScale**2 - sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)/bScale**2)*sign(Cartz)),
   *    xx[2] = atan2(Carty, Cartx)]"
   */
  const double tmp_1 = ((bScale)*(bScale));
  const double tmp_2 = ((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz));
  const double tmp_3 = sqrt(-4*((Cartz)*(Cartz))*tmp_1 + ((bScale)*(bScale)*(bScale)*(bScale)) + 2*tmp_1*tmp_2 + ((tmp_2)*(tmp_2)));
  const double tmp_4 = (1.0/(tmp_1));
  xx[0] = M_SQRT1_2*sqrt(-tmp_1 + tmp_2 + tmp_3);
  xx[1] = acos(M_SQRT1_2*sqrt(tmp_2*tmp_4 - tmp_3*tmp_4 + 1)*(((Cartz) > 0) - ((Cartz) < 0)));
  xx[2] = atan2(Carty, Cartx);

    // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
    // Recall that:
    // xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0; // Cell-centered grid.
    //   --> j = (int) ( (xx[0][j] - xxmin[0]) / params->dxx0 + (1.0/2.0) + NGHOSTS )
    Cart_to_i0i1i2[0] = (int)( ( xx[0] - (0) ) / params->dxx0 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
    Cart_to_i0i1i2[1] = (int)( ( xx[1] - (0) ) / params->dxx1 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
    Cart_to_i0i1i2[2] = (int)( ( xx[2] - (-M_PI) ) / params->dxx2 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
}
