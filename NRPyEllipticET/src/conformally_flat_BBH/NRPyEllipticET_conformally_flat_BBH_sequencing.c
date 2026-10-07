#include "../NRPyEllipticET.h"
#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"

/*
 * Grid sequencing for the leapfrog hyperbolic relaxation (EXPERIMENTAL; coarse_grid_levels > 0).
 *
 * Relax first on grids with Nxx/2 (and Nxx/4) points per direction: each step there is ~8x cheaper
 * and the time step ~2x larger. Then interpolate to the next finer grid and continue. The leapfrog
 * relaxation damps all modes at roughly the same rate, so interpolation error (whose residual is
 * nearly as large as the initial one) would otherwise cost almost a full solve; hence, after each
 * interpolation, the first prolongation_smoothing_steps steps use eta = 2/dt, for which the leapfrog
 * update reduces to the Richardson iteration uu^{n+1} = uu^n + (dt^2/2) c^2 R(uu^n), which removes
 * it quickly (like a multigrid smoother). The final answer is the 10th-order solution on the full grid.
 *
 * The 4th-order-first variant (cheaper stencils at the start) was tried and dropped: a 4th-order step
 * costs about as much as a 10th-order one, since most of the RHS cost is not in the stencils.
 */

// Set up a complete griddata on an Nxx grid. Physical parameters come from params_template;
// the grid-dependent ones are recomputed. Mirrors the setup in the relaxation driver.
void conformally_flat_BBH_sequencing_griddata_setup(const paramstruct *restrict params_template, const int Nxx[3],
                                                    const REAL CFL_FACTOR, const REAL eta, griddata_struct *restrict gd,
                                                    REAL *restrict dt) {
  gd->params = *params_template;
  {
    const int EigenCoord = 1;
    conformally_flat_BBH_set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &gd->params, gd->xx);
    conformally_flat_BBH_driver_bcstruct(&gd->params, &gd->bcstruct, gd->xx);
    for (int i = 0; i < 3; i++) free(gd->xx[i]);
  }
  conformally_flat_BBH_set_Nxx_dxx_invdx_params__and__xx(0, Nxx, &gd->params, gd->xx);
  *dt = conformally_flat_BBH_find_timestep(&gd->params, gd->xx, CFL_FACTOR);
  conformally_flat_BBH_MoL_malloc_y_n_gfs(&gd->params, &gd->gridfuncs);
  conformally_flat_BBH_rfm_precompute_rfmstruct_malloc(&gd->params, &gd->rfmstruct);
  conformally_flat_BBH_rfm_precompute_rfmstruct_define(&gd->params, gd->xx, &gd->rfmstruct);
  conformally_flat_BBH_initial_guess_all_points(&gd->params, gd->xx, gd->gridfuncs.y_n_gfs);
  conformally_flat_BBH_MoL_malloc_non_y_n_gfs(&gd->params, &gd->gridfuncs);
  conformally_flat_BBH_auxevol_gfs_all_points(&gd->params, gd->xx, gd->gridfuncs.auxevol_gfs);
  conformally_flat_BBH_wavespeed_gf_all_points(&gd->params, CFL_FACTOR, *dt, gd->xx, gd->gridfuncs.auxevol_gfs);
  gd->params.eta_damping = eta;
}

void conformally_flat_BBH_sequencing_griddata_free(griddata_struct *restrict gd) {
  conformally_flat_BBH_rfm_precompute_rfmstruct_freemem(&gd->params, &gd->rfmstruct);
  conformally_flat_BBH_freemem_bcstruct(&gd->params, &gd->bcstruct);
  conformally_flat_BBH_MoL_free_memory_y_n_gfs(&gd->params, &gd->gridfuncs);
  conformally_flat_BBH_MoL_free_memory_non_y_n_gfs(&gd->params, &gd->gridfuncs);
  for (int i = 0; i < 3; i++) free(gd->xx[i]);
}

// Leapfrog relaxation of gd from rest (uu^{-1} = uu^0) until log10(L2 residual within residual_radius)
// < log_target, checked every info_output_freq steps (every step if 0), like the standard loop. The
// first n_smooth steps use eta = 2/dt (Richardson smoothing; see above). Returns the number of steps,
// or -1 if max_its is reached or the residual becomes non-finite.
int conformally_flat_BBH_sequencing_relax(griddata_struct *restrict gd, const REAL dt, const int n_smooth, const REAL log_target,
                                          const int max_its, const REAL residual_radius, REAL *restrict final_log_residual) {
  DECLARE_CCTK_PARAMETERS;
  const int Ntot = gd->params.Nxx_plus_2NGHOSTS0 * gd->params.Nxx_plus_2NGHOSTS1 * gd->params.Nxx_plus_2NGHOSTS2;
  REAL *restrict weights = conformally_flat_BBH_L2_norm_weights_malloc(&gd->params, residual_radius, gd->xx);
  if (weights == NULL) CCTK_ERROR("Could not allocate L2-norm weights.");
  memcpy(gd->gridfuncs.k_even_gfs, gd->gridfuncs.y_n_gfs, sizeof(REAL) * NUM_EVOL_GFS * Ntot);
  const REAL eta = gd->params.eta_damping;
  const int check_every = info_output_freq > 0 ? info_output_freq : 1;
  int result = -1;
  *final_log_residual = NAN;
  for (int n = 0; n < max_its; n++) {
    gd->params.eta_damping = (n < n_smooth) ? 2.0 / dt : eta;
    conformally_flat_BBH_leapfrog_step_forward_in_time(gd, dt);
    gd->params.eta_damping = eta;
    if (n % check_every != 0) continue;
    // Residual of the time level before this step (computed by the leapfrog step).
    const REAL log_res = conformally_flat_BBH_L2_norm_of_gf_weighted(&gd->params, UUGF, weights, gd->gridfuncs.diagnostic_output_gfs);
    *final_log_residual = log_res;
    if (verbose)
      CCTK_VINFO("  relax %dx%dx%d it %6d log10 residual %.3f", gd->params.Nxx0, gd->params.Nxx1, gd->params.Nxx2, n, log_res);
    if (!isfinite(log_res)) break;
    if (log_res < log_target) {
      result = n + 1;
      break;
    }
  }
  free(weights);
  return result;
}

// Cubic Lagrange weights for nodes 0,1,2,3 evaluated at t.
static void cubic_weights(const REAL t, REAL w[4]) {
  w[0] = -(t - 1.0) * (t - 2.0) * (t - 3.0) / 6.0;
  w[1] = t * (t - 2.0) * (t - 3.0) / 2.0;
  w[2] = -t * (t - 1.0) * (t - 3.0) / 2.0;
  w[3] = t * (t - 1.0) * (t - 2.0) / 6.0;
}

// Base index and weights of a 4-point stencil around coordinate x on a uniform coordinate array.
static void stencil_1d(const REAL *restrict xc, const int Nc, const REAL dxc, const REAL x, int *restrict base, REAL w[4]) {
  const REAL p = (x - xc[0]) / dxc;
  int b = (int)floor(p) - 1;
  if (b < 0) b = 0;
  if (b > Nc - 4) b = Nc - 4;
  *base = b;
  cubic_weights(p - (REAL)b, w);
}

// Interpolate uu from the coarse grid c to every point (including ghost zones) of the fine grid f,
// by tensor-product cubic interpolation in the uniform reference coordinates (xx0, xx1, xx2).
// The coarse ghost zones (valid through the relaxation's boundary treatment) cover the fine ones.
void conformally_flat_BBH_sequencing_prolongate(const griddata_struct *restrict c, griddata_struct *restrict f) {
  const int Nc0 = c->params.Nxx_plus_2NGHOSTS0, Nc1 = c->params.Nxx_plus_2NGHOSTS1, Nc2 = c->params.Nxx_plus_2NGHOSTS2;
  const int Nf0 = f->params.Nxx_plus_2NGHOSTS0, Nf1 = f->params.Nxx_plus_2NGHOSTS1, Nf2 = f->params.Nxx_plus_2NGHOSTS2;
  const REAL *restrict uc = c->gridfuncs.y_n_gfs; // UUGF = 0: first Nc0*Nc1*Nc2 entries
  REAL *restrict uf = f->gridfuncs.y_n_gfs;
#pragma omp parallel for
  for (int i2 = 0; i2 < Nf2; i2++) {
    int b2; REAL w2[4];
    stencil_1d(c->xx[2], Nc2, c->params.dxx2, f->xx[2][i2], &b2, w2);
    for (int i1 = 0; i1 < Nf1; i1++) {
      int b1; REAL w1[4];
      stencil_1d(c->xx[1], Nc1, c->params.dxx1, f->xx[1][i1], &b1, w1);
      for (int i0 = 0; i0 < Nf0; i0++) {
        int b0; REAL w0[4];
        stencil_1d(c->xx[0], Nc0, c->params.dxx0, f->xx[0][i0], &b0, w0);
        REAL sum = 0.0;
        for (int k = 0; k < 4; k++)
          for (int j = 0; j < 4; j++)
            for (int i = 0; i < 4; i++)
              sum += w2[k] * w1[j] * w0[i] * uc[(b0 + i) + Nc0 * ((b1 + j) + Nc1 * (b2 + k))];
        uf[i0 + Nf0 * (i1 + Nf1 * i2)] = sum;
      } // END LOOP: i0
    } // END LOOP: i1
  } // END LOOP: i2
}
