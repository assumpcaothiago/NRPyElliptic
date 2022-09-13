void conformally_flat_BBH_MoL_malloc_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void conformally_flat_BBH_MoL_free_memory_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void conformally_flat_BBH_MoL_malloc_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void conformally_flat_BBH_MoL_free_memory_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void conformally_flat_BBH_MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt);
void conformally_flat_BBH_set_bcstruct(const paramstruct *restrict params,
                                       gz_map *restrict bc_gz_map,
                                       parity_condition *bc_parity_conditions,
                                       bc_struct *restrict bcstruct);
void conformally_flat_BBH_driver_bcstruct(const paramstruct *restrict params, bc_struct *restrict bcstruct, REAL *restrict xx[3]);
void conformally_flat_BBH_freemem_bcstruct(const paramstruct *restrict params, const bc_struct *restrict bcstruct);
void conformally_flat_BBH_set_up__bc_gz_map_and_parity_condns(const paramstruct *restrict params,
                                                              REAL *restrict xx[3], gz_map *restrict bc_gz_map,parity_condition *restrict bc_parity_conditions);
REAL conformally_flat_BBH_find_timestep(const paramstruct *restrict params, REAL *restrict xx[3], const REAL CFL_FACTOR);
void conformally_flat_BBH_xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]);
void conformally_flat_BBH_set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],paramstruct *restrict params, REAL *restrict xx[3]);
void conformally_flat_BBH_rfm_precompute_rfmstruct_malloc(const paramstruct *restrict params, rfm_struct *restrict rfmstruct);
void conformally_flat_BBH_rfm_precompute_rfmstruct_define(const paramstruct *restrict params, REAL *restrict xx[3], rfm_struct *restrict rfmstruct);
void conformally_flat_BBH_rfm_precompute_rfmstruct_freemem(const paramstruct *restrict params, rfm_struct *restrict rfmstruct);
void conformally_flat_BBH_set_Cparameters_to_default(paramstruct *restrict params);
void conformally_flat_BBH_initial_guess_single_point(const paramstruct *restrict params,
                                                     const REAL xx0, const REAL xx1, const REAL xx2,
                                                     REAL *uu_exact, REAL *vv_exact);
void conformally_flat_BBH_initial_guess_all_points(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs);
void conformally_flat_BBH_rhs_eval(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs);
void conformally_flat_BBH_auxevol_gfs_all_points(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs);
void conformally_flat_BBH_wavespeed_gf_all_points(const paramstruct *restrict params, const REAL CFL_FACTOR, const REAL dt, REAL *restrict xx[3], REAL *restrict auxevol_gfs);
void conformally_flat_BBH_residual_all_points(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict residual_gf);
REAL conformally_flat_BBH_L2_norm_of_gf(const paramstruct *restrict params, const int gf_index, const REAL integration_radius, REAL *restrict xx[3], const REAL *restrict in_gf);
void conformally_flat_BBH_gridfunction_z_axis(const paramstruct *restrict params, const int gf_index, REAL *restrict xx[3], const REAL *restrict in_gf, FILE *restrict outfile);
void conformally_flat_BBH_compute_ADM_Cartesian_quantities_from_uu(const cGH *restrict cctkGH,
                                                                   const char *restrict orbital_plane,
                                                                   const paramstruct *restrict params,
                                                                   const REAL *restrict x,
                                                                   const REAL *restrict y,
                                                                   const REAL *restrict z,
                                                                   const REAL *restrict position_shift,
                                                                   const REAL *restrict in_gfs,
                                                                   REAL *restrict alpha,
                                                                   REAL *restrict betaU0, REAL *restrict betaU1, REAL *restrict betaU2,
                                                                   REAL *restrict gammaDD00, REAL *restrict gammaDD01, REAL *restrict gammaDD02,
                                                                   REAL *restrict gammaDD11, REAL *restrict gammaDD12,
                                                                   REAL *restrict gammaDD22,
                                                                   REAL *restrict KDD00, REAL *restrict KDD01, REAL *restrict KDD02,
                                                                   REAL *restrict KDD11, REAL *restrict KDD12,
                                                                   REAL *restrict KDD22);
void conformally_flat_BBH_get_xx_from_xyz(const double AMAX, const double SINHWAA, const double bScale,
                                          const double xCart[3], double xx[3]);
void conformally_flat_BBH_lapse_averaged_initial_data(const cGH *restrict cctkGH,
                                                      const char *restrict orbital_plane,
                                                      const paramstruct *restrict params,
                                                      const REAL *restrict x,
                                                      const REAL *restrict y,
                                                      const REAL *restrict z,
                                                      const REAL *restrict position_shift,
                                                      const REAL *restrict in_gfs,
                                                      REAL *restrict alpha);
void conformally_flat_BBH_lapse_power_of_psi_initial_data(const cGH *restrict cctkGH,
                                                          const paramstruct *restrict params,
                                                          const REAL *restrict x,
                                                          const REAL *restrict y,
                                                          const REAL *restrict z,
                                                          const REAL *restrict in_gfs,
                                                          REAL *restrict alpha);
void conformally_flat_BBH_set_optimum_eta_damping(const CCTK_REAL dt,
                                                  const CCTK_REAL CFL_FACTOR,
                                                  paramstruct *restrict params);

void NRPyEllipticET_conformally_flat_BBH_interpolate_solution_to_ADMBase( const CCTK_INT num_pts,
                                                                          const CCTK_INT *input_array_dims,
                                                                          const CCTK_REAL *origin,
                                                                          const CCTK_REAL *deltas,
                                                                          const void **interp_coords,
                                                                          const CCTK_REAL *input_gf,
                                                                          CCTK_REAL *output_gf );
