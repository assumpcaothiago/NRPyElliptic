#ifndef TEST_CCTK_PARAMETERS_H
#define TEST_CCTK_PARAMETERS_H
#define DECLARE_CCTK_PARAMETERS
extern int N0, N1, N2, max_iterations, info_output_freq, verbose;
extern int coarse_grid_levels, prolongation_smoothing_steps;
extern double CFL_FACTOR, eta_damping, log_target_residual;
extern double coarse_log_target_residual, residual_integration_radius;
extern double domain_size, sinh_width, foci_position, number_of_LCT;
extern double conformally_flat_BBH_puncture_0_bare_mass, conformally_flat_BBH_puncture_1_bare_mass;
extern double conformally_flat_BBH_puncture_0_pos[3], conformally_flat_BBH_puncture_1_pos[3];
extern double conformally_flat_BBH_puncture_0_P[3], conformally_flat_BBH_puncture_1_P[3];
extern double conformally_flat_BBH_puncture_0_S[3], conformally_flat_BBH_puncture_1_S[3];
extern double lapse_exponent_n;
extern const char *orbital_plane, *time_integrator, *initial_lapse;
#endif
