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
  printf("puncture_0_bare_mass  = %lf\n", puncture_0_bare_mass);
  printf("puncture_1_bare_mass  = %lf\n", puncture_1_bare_mass);
  printf("puncture_0_x = %lf\n", puncture_0_x);
  printf("puncture_0_y = %lf\n", puncture_0_y);
  printf("puncture_0_z = %lf\n", puncture_0_z);
  printf("puncture_1_x = %lf\n", puncture_1_x);
  printf("puncture_1_y = %lf\n", puncture_1_y);
  printf("puncture_1_z = %lf\n", puncture_1_z);
  printf("puncture_0_P_x         = %lf\n", puncture_0_P_x);
  printf("puncture_0_P_y         = %lf\n", puncture_0_P_y);
  printf("puncture_0_P_z         = %lf\n", puncture_0_P_z);
  printf("puncture_1_P_x         = %lf\n", puncture_1_P_x);
  printf("puncture_1_P_y         = %lf\n", puncture_1_P_y);
  printf("puncture_1_P_z         = %lf\n", puncture_1_P_z);
  printf("puncture_0_S_x         = %lf\n", puncture_0_S_x);
  printf("puncture_0_S_y         = %lf\n", puncture_0_S_y);
  printf("puncture_0_S_z         = %lf\n", puncture_0_S_z);
  printf("puncture_1_S_x         = %lf\n", puncture_1_S_x);
  printf("puncture_1_S_y         = %lf\n", puncture_1_S_y);
  printf("puncture_1_S_z         = %lf\n", puncture_1_S_z);
  printf("Nxx0         = %d\n",  Nxx0);
  printf("Nxx1         = %d\n",  Nxx1);
  printf("Nxx2         = %d\n",  Nxx2);
  printf("NGHOSTS      = %d\n",  NGHOSTS);
}
