void MoL_malloc_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_malloc_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt);
void set_bcstruct(const paramstruct *restrict params,
                  gz_map *restrict bc_gz_map,
                  parity_condition *bc_parity_conditions,
                  bc_struct *restrict bcstruct);
void driver_bcstruct(const paramstruct *restrict params, bc_struct *restrict bcstruct, REAL *restrict xx[3]);
void apply_bcs_curvilinear_inner_only(const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict xx[3],
                           REAL *restrict gfs);
void freemem_bcstruct(const paramstruct *restrict params, const bc_struct *restrict bcstruct);
void set_up__bc_gz_map_and_parity_condns(const paramstruct *restrict params,
                                             REAL *restrict xx[3], gz_map *restrict bc_gz_map,parity_condition *restrict bc_parity_conditions);
void apply_bcs_curvilinear_radiation(const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict xx[3],
                           REAL *restrict gfs, REAL *restrict gfs_rhss);
REAL find_timestep(const paramstruct *restrict params, REAL *restrict xx[3], const REAL CFL_FACTOR);
void xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]);
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],paramstruct *restrict params, REAL *restrict xx[3]);
void Cart_to_xx_and_nearest_i0i1i2(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]);
void Cart_to_xx_and_nearest_i0i1i2_global_grid_center(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]);
void rfm_precompute_rfmstruct_malloc(const paramstruct *restrict params, rfm_struct *restrict rfmstruct);
void rfm_precompute_rfmstruct_define(const paramstruct *restrict params, REAL *restrict xx[3], rfm_struct *restrict rfmstruct);
void rfm_precompute_rfmstruct_freemem(const paramstruct *restrict params, rfm_struct *restrict rfmstruct);
void set_Cparameters_to_default(paramstruct *restrict params);
void initial_guess_single_point(const paramstruct *restrict params,
                const REAL xx0, const REAL xx1, const REAL xx2,
                REAL *uu_exact, REAL *vv_exact);
void initial_guess_all_points(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs);
void rhs_eval(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs);
void auxevol_gfs_all_points(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs);
void wavespeed_gf_all_points(const paramstruct *restrict params, const REAL CFL_FACTOR, const REAL dt, REAL *restrict xx[3], REAL *restrict auxevol_gfs);
void residual_all_points(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict residual_gf);
REAL L2_norm_of_gf(const paramstruct *restrict params, const int gf_index, const REAL integration_radius, REAL *restrict xx[3], const REAL *restrict in_gf);
void gridfunction_z_axis(const paramstruct *restrict params, const int gf_index, REAL *restrict xx[3], const REAL *restrict in_gf, FILE *restrict outfile);
void print_puncture_parameters(const paramstruct *restrict params);
int main(int argc, const char *argv[]);
void compute_ADM_Cartesian_quantities_from_uu(const cGH *restrict cctkGH,
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
void get_xx_from_xyz(const double AMAX, const double SINHWAA, const double bScale,
                     const double xCart[3], double xx[3]);
void lapse_averaged_initial_data(const cGH *restrict cctkGH,
                                 const char *restrict orbital_plane,
                                 const paramstruct *restrict params,
                                 const REAL *restrict x,
                                 const REAL *restrict y,
                                 const REAL *restrict z,
                                 const REAL *restrict position_shift,
                                 const REAL *restrict in_gfs,
                                 REAL *restrict alpha);
void lapse_power_of_psi_initial_data(const cGH *restrict cctkGH,
                                     const paramstruct *restrict params,
                                     const REAL *restrict x,
                                     const REAL *restrict y,
                                     const REAL *restrict z,
                                     const REAL *restrict in_gfs,
                                     REAL *restrict alpha);
void NRPyEllipticET_cleanup(CCTK_ARGUMENTS);
void NRPyEllipticET_interpolate_solution_to_ADMBase( const CCTK_INT num_pts,
                                                     const CCTK_INT input_array_dims[3],
                                                     const CCTK_REAL origin[3],
                                                     const CCTK_REAL deltas[3],
                                                     const void *restrict interp_coords[3],
                                                     const CCTK_REAL *restrict input_gf,
                                                     CCTK_REAL *restrict output_gf );
void NRPyEllipticET_PunctureInitialData_Hyperbolic_Relaxation();
void NRPyEllipticET_PunctureInitialData_Initialize_ADMBase( const cGH *restrict cctkGH,
                                                            const int *restrict cctk_lsh,
                                                            const CCTK_REAL *restrict x,
                                                            const CCTK_REAL *restrict y,
                                                            const CCTK_REAL *restrict z,
                                                            CCTK_REAL *restrict uuGF,
                                                            CCTK_REAL *restrict alp,
                                                            CCTK_REAL *restrict betax,
                                                            CCTK_REAL *restrict betay,
                                                            CCTK_REAL *restrict betaz,
                                                            CCTK_REAL *restrict gxx,
                                                            CCTK_REAL *restrict gxy,
                                                            CCTK_REAL *restrict gxz,
                                                            CCTK_REAL *restrict gyy,
                                                            CCTK_REAL *restrict gyz,
                                                            CCTK_REAL *restrict gzz,
                                                            CCTK_REAL *restrict kxx,
                                                            CCTK_REAL *restrict kxy,
                                                            CCTK_REAL *restrict kxz,
                                                            CCTK_REAL *restrict kyy,
                                                            CCTK_REAL *restrict kyz,
                                                            CCTK_REAL *restrict kzz );
void NRPyEllipticET(CCTK_ARGUMENTS);
void NRPyEllipticET_init(CCTK_ARGUMENTS);
