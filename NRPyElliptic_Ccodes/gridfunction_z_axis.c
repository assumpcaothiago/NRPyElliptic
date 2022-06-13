#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Output gridfunction along z axis
 */
void gridfunction_z_axis(const paramstruct *restrict params, const int gf_index, REAL *restrict xx[3], const REAL *restrict in_gf, FILE *restrict outfile) {
#include "./set_Cparameters.h"


  // Set i2 to any index (z-axis is independent of i2)
  const int i2 = Nxx_plus_2NGHOSTS2/2; const REAL xx2 = xx[2][i2];
  
  // Output postive z-axis
  // Set x1 = 0 (i1 = NGHOSTS)
  int i1 = NGHOSTS; REAL xx1 = xx[1][i1];
  for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++){
    REAL xx0 = xx[0][i0];
    const double tmp_0 = (1.0/(SINHWAA));
    const REAL z = sqrt(((AMAX)*(AMAX))*((exp(tmp_0*xx0) - exp(-tmp_0*xx0))*(exp(tmp_0*xx0) - exp(-tmp_0*xx0)))/((exp(tmp_0) - exp(-tmp_0))*(exp(tmp_0) - exp(-tmp_0))) + ((bScale)*(bScale)))*cos(xx1);

    const REAL gf_at_z = in_gf[IDX4S(gf_index, i0, i1, i2)];
    fprintf(outfile, "%.17e %.17e\n", z, gf_at_z);

  } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)

  // Output negative z-axis
  // Set x1 = pi (i1 = Nxx_plus_2NGHOSTS1 - NGHOSTS)
  i1 = Nxx_plus_2NGHOSTS1 - NGHOSTS - 1; xx1 = xx[1][i1];
  for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++){
    REAL xx0 = xx[0][i0];
    const double tmp_0 = (1.0/(SINHWAA));
    const REAL z = sqrt(((AMAX)*(AMAX))*((exp(tmp_0*xx0) - exp(-tmp_0*xx0))*(exp(tmp_0*xx0) - exp(-tmp_0*xx0)))/((exp(tmp_0) - exp(-tmp_0))*(exp(tmp_0) - exp(-tmp_0))) + ((bScale)*(bScale)))*cos(xx1);

    const REAL gf_at_z = in_gf[IDX4S(gf_index, i0, i1, i2)];
    fprintf(outfile, "%.17e %.17e\n", z, gf_at_z);

  } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)

  // Add i0 = NGHOSTS line (x0 = 0)
  // Here, we do not add the points i1 = NGHOSTS and i1 = xx_plus_2NGHOSTS1 - NGHOSTS,
  // as they hare already been added in the previous loops
  int i0 = NGHOSTS; REAL xx0 = xx[0][i0];
  for (int i1 = NGHOSTS + 1; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS - 1; i1++){
    REAL xx1 = xx[1][i1];
    const double tmp_0 = (1.0/(SINHWAA));
    const REAL z = sqrt(((AMAX)*(AMAX))*((exp(tmp_0*xx0) - exp(-tmp_0*xx0))*(exp(tmp_0*xx0) - exp(-tmp_0*xx0)))/((exp(tmp_0) - exp(-tmp_0))*(exp(tmp_0) - exp(-tmp_0))) + ((bScale)*(bScale)))*cos(xx1);

    const REAL gf_at_z = in_gf[IDX4S(gf_index, i0, i1, i2)];
    fprintf(outfile, "%.17e %.17e\n", z, gf_at_z);

  } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
}
