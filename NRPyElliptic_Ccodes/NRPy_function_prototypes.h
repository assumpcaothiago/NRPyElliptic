void MoL_malloc_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_malloc_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt);
void bcstruct_set_up(const paramstruct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct);
void apply_bcs_outerradiation_and_inner(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                        const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS],
                                        REAL *restrict gfs, REAL *restrict rhs_gfs);
void apply_bcs_inner_only(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
void apply_bcs_outerextrap_and_inner(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
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
void gridfunction_xz_plane(const paramstruct *restrict params, const int gf_index, REAL *restrict xx[3], const REAL *restrict in_gf, FILE *restrict outfile);
void print_puncture_parameters(const paramstruct *restrict params);
REAL compute_wavespeed_at_OB(const paramstruct *restrict params, REAL *restrict auxevol_gfs);
int main(int argc, const char *argv[]);
