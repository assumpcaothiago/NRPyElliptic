#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Print all puncture parameters
 */
void print_puncture_parameters(const paramstruct *restrict params) {
#include "./set_Cparameters.h"


  printf("bScale       = %lf\n", bScale);
  printf("SINHWAA      = %lf\n", SINHWAA);
  printf("AMAX         = %lf\n", AMAX);
  printf("eta_damping  = %lf\n", eta_damping);
  printf("bare_mass_0  = %lf\n", bare_mass_0);
  printf("bare_mass_1  = %lf\n", bare_mass_1);
  printf("puncture_0_x = %lf\n", puncture_0_x);
  printf("puncture_0_y = %lf\n", puncture_0_y);
  printf("puncture_0_z = %lf\n", puncture_0_z);
  printf("puncture_1_x = %lf\n", puncture_1_x);
  printf("puncture_1_y = %lf\n", puncture_1_y);
  printf("puncture_1_z = %lf\n", puncture_1_z);
  printf("P0_x         = %lf\n", P0_x);
  printf("P0_y         = %lf\n", P0_y);
  printf("P0_z         = %lf\n", P0_z);
  printf("P1_x         = %lf\n", P1_x);
  printf("P1_y         = %lf\n", P1_y);
  printf("P1_z         = %lf\n", P1_z);
  printf("S0_x         = %lf\n", S0_x);
  printf("S0_y         = %lf\n", S0_y);
  printf("S0_z         = %lf\n", S0_z);
  printf("S1_x         = %lf\n", S1_x);
  printf("S1_y         = %lf\n", S1_y);
  printf("S1_z         = %lf\n", S1_z);
  printf("Nxx0         = %d\n",  Nxx0);
  printf("Nxx1         = %d\n",  Nxx1);
  printf("Nxx2         = %d\n",  Nxx2);
  printf("NGHOSTS      = %d\n",  NGHOSTS);
}
