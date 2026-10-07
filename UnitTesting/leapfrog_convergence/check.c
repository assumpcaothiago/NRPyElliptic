/* Exercise the real C relaxation drivers without an Einstein Toolkit install. */
#include "NRPyEllipticET.h"
#include "conformally_flat_BBH_NRPy_basic_defines.h"
#include "conformally_flat_BBH_NRPy_function_prototypes.h"

int N0 = 48, N1 = 48, N2 = 12, max_iterations = 10000, info_output_freq = 0, verbose = 0;
int coarse_grid_levels = 0, prolongation_smoothing_steps = 400;
double CFL_FACTOR = 0.5, eta_damping = 6, log_target_residual = -5;
double coarse_log_target_residual = -6, residual_integration_radius = 100;
double domain_size = 1e6, sinh_width = 0.07, foci_position = 5, number_of_LCT = 1;
double conformally_flat_BBH_puncture_0_bare_mass = 0.51841937253070269;
double conformally_flat_BBH_puncture_1_bare_mass = 0.39192387717971311;
double conformally_flat_BBH_puncture_0_pos[3] = {0, 0, 5};
double conformally_flat_BBH_puncture_1_pos[3] = {0, 0, -5};
double conformally_flat_BBH_puncture_0_P[3] = {0.09530152296974252, 0, -0.00084541526517121};
double conformally_flat_BBH_puncture_1_P[3] = {-0.09530152296974252, 0, 0.00084541526517121};
double conformally_flat_BBH_puncture_0_S[3] = {0, 0.09509112426035504, 0};
double conformally_flat_BBH_puncture_1_S[3] = {0, -0.09156449704142013, 0};
double lapse_exponent_n = -2;
const char *orbital_plane = "xz", *time_integrator = "leapfrog", *initial_lapse = "psi^n";
CCTK_REAL *NRPyEllipticET_uu, *NRPyEllipticET_xx[3];

static paramstruct parameters(void) {
  griddata_struct griddata;
  conformally_flat_BBH_set_Cparameters_to_default(&griddata.params);
#include "conformally_flat_BBH_set_Cparameters_from_parfile.h"
  griddata.params.eta_damping = eta_damping;
  return griddata.params;
}

/* Recompute from the returned field and use the original, uncached norm. */
static double actual_residual(griddata_struct *gd) {
  conformally_flat_BBH_residual_all_points(&gd->params, &gd->rfmstruct, gd->gridfuncs.auxevol_gfs,
                                          gd->gridfuncs.y_n_gfs, gd->gridfuncs.diagnostic_output_gfs);
  return conformally_flat_BBH_L2_norm_of_gf(&gd->params, UUGF, residual_integration_radius,
                                           gd->xx, gd->gridfuncs.diagnostic_output_gfs);
}

static void free_driver_output(void) {
  free(NRPyEllipticET_uu);
  for (int i = 0; i < 3; i++) free(NRPyEllipticET_xx[i]);
}

int main(int argc, char **argv) {
  if (argc != 2) return 1;
  const int final_check = strcmp(argv[1], "sequencing-final") == 0;
  const int standard = strcmp(argv[1], "standard-cached") == 0;
  if (strcmp(argv[1], "standard-exhausted") == 0) {
    max_iterations = 1;
    log_target_residual = -14;
    /* This must raise the production driver's fatal nonconvergence error. */
    NRPyEllipticET_conformally_flat_BBH_Hyperbolic_Relaxation();
    free_driver_output();
    return 0;
  }
  if (!final_check && !standard && strcmp(argv[1], "sequencing-cached") != 0) return 1;
  if (final_check) {
    eta_damping = 12.5;
    log_target_residual = -12;
    info_output_freq = 1000;
    max_iterations = 2000;
  }
  const paramstruct params = parameters();
  const int Nxx[3] = {N0, N1, N2};
  griddata_struct gd;
  double dt, reported = NAN;
  int steps = 1;
  if (standard) NRPyEllipticET_conformally_flat_BBH_Hyperbolic_Relaxation();
  conformally_flat_BBH_sequencing_griddata_setup(&params, Nxx, CFL_FACTOR, eta_damping, &gd, &dt);
  if (standard) {
    const int count = gd.params.Nxx_plus_2NGHOSTS0 * gd.params.Nxx_plus_2NGHOSTS1 * gd.params.Nxx_plus_2NGHOSTS2;
    memcpy(gd.gridfuncs.y_n_gfs, NRPyEllipticET_uu, count * sizeof(REAL));
    free_driver_output();
  } else {
    if (final_check) {
      const int Nc[3] = {N0 / 2, N1 / 2, N2 / 2};
      griddata_struct coarse;
      double dtc, coarse_residual;
      conformally_flat_BBH_sequencing_griddata_setup(&params, Nc, CFL_FACTOR, eta_damping, &coarse, &dtc);
      const int coarse_steps = conformally_flat_BBH_sequencing_relax(&coarse, dtc, 0,
          coarse_log_target_residual, 10000, residual_integration_radius, &coarse_residual);
      if (coarse_steps < 0) {
        fprintf(stderr, "Coarse-grid setup did not converge\n");
        conformally_flat_BBH_sequencing_griddata_free(&coarse);
        conformally_flat_BBH_sequencing_griddata_free(&gd);
        return 1;
      }
      conformally_flat_BBH_sequencing_prolongate(&coarse, &gd);
      conformally_flat_BBH_sequencing_griddata_free(&coarse);
    }
    steps = conformally_flat_BBH_sequencing_relax(&gd, dt,
        final_check ? prolongation_smoothing_steps : 0, log_target_residual,
        max_iterations, residual_integration_radius, &reported);
  }
  const double actual = actual_residual(&gd);
  const int passed = steps > 0 && isfinite(actual) && actual < log_target_residual &&
      (standard || (isfinite(reported) && fabs(actual - reported) < 1e-10)) &&
      (!final_check || steps == max_iterations);
  printf("%s: target=%.2f actual=%.16g", argv[1], log_target_residual, actual);
  if (!standard) printf(" reported=%.16g steps=%d", reported, steps);
  printf("\n");
  conformally_flat_BBH_sequencing_griddata_free(&gd);
  return passed ? 0 : 1;
}
