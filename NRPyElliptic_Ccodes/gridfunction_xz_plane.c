#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Output gridfunction along xz-plane
 */
void gridfunction_xz_plane(const paramstruct *restrict params, const int gf_index, REAL *restrict xx[3], const REAL *restrict in_gf, FILE *restrict outfile) {
#include "./set_Cparameters.h"


  // Set x2 = -pi (i2 = NGHOSTS) and loop over i0 and i1
  int i2 = NGHOSTS; REAL xx2 = xx[2][i2];
  
  // Output points with x < 0
  for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++){
    REAL xx1 = xx[1][i1];
    for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++){
      REAL xx0 = xx[0][i0];
      const double tmp_0 = (1.0/(SINHWAA));
      const double tmp_1 = exp(tmp_0) - exp(-tmp_0);
      const double tmp_3 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
      const REAL xCart = AMAX*tmp_3*sin(xx1)*cos(xx2)/tmp_1;
      const REAL zCart = sqrt(((AMAX)*(AMAX))*((tmp_3)*(tmp_3))/((tmp_1)*(tmp_1)) + ((bScale)*(bScale)))*cos(xx1);

      const REAL gf_at_xz = in_gf[IDX4S(gf_index, i0, i1, i2)];
      fprintf(outfile, "%.17e %.17e %.17e\n", xCart, zCart, gf_at_xz);

    } // END for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
  } // END for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)

  // Set x2 = 0 (i2 = Nxx_plus_2NGHOSTS2/2) and loop over i0 and i1
  i2 = Nxx_plus_2NGHOSTS2/2; xx2 = xx[2][i2];
  
  // Output points with x > 0
  for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++){
    REAL xx1 = xx[1][i1];
    for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++){
      REAL xx0 = xx[0][i0];
      const double tmp_0 = (1.0/(SINHWAA));
      const double tmp_1 = exp(tmp_0) - exp(-tmp_0);
      const double tmp_3 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
      const REAL xCart = AMAX*tmp_3*sin(xx1)*cos(xx2)/tmp_1;
      const REAL zCart = sqrt(((AMAX)*(AMAX))*((tmp_3)*(tmp_3))/((tmp_1)*(tmp_1)) + ((bScale)*(bScale)))*cos(xx1);

      const REAL gf_at_xz = in_gf[IDX4S(gf_index, i0, i1, i2)];
      fprintf(outfile, "%.17e %.17e %.17e\n", xCart, zCart, gf_at_xz);

    } // END for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
  } // END for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
}
