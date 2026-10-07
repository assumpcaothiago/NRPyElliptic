#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Leapfrog alternative to conformally_flat_BBH_MoL_step_forward_in_time().
 *
 * Eliminating vv from the relaxation system uu_t = vv - eta*uu, vv_t = c^2 R(uu)
 * gives the damped wave equation uu_tt + eta*uu_t = c^2 R(uu), where R is the
 * residual of the elliptic equation. Discretizing both time derivatives with
 * centered differences gives
 *
 *   uu^{n+1} = ( 2 uu^n - (1 - eta dt/2) uu^{n-1} + dt^2 c^2 R(uu^n) ) / (1 + eta dt/2),
 *
 * which needs one RHS evaluation per step (RK4 needs four), and whose steady
 * state is identical to RK4's (all time-derivative terms vanish). It is stable
 * for local Courant factors up to ~0.71 of the RK4 limit, so the effective
 * Courant factor (Courant factor / RHS evaluations per step) is ~2.8x larger.
 *
 * Storage: y_n_gfs[UUGF] holds uu^n and k_even_gfs[UUGF] holds uu^{n-1}; the two
 * pointers are swapped on exit, so y_n_gfs[UUGF] always holds the newest uu.
 * k_odd_gfs is scratch. Since rhs_eval already computes c^2 R(uu^n), the residual
 * of uu^n is stored in diagnostic_output_gfs[UUGF] at no extra cost.
 *
 * Outer ghost zones: the Sommerfeld condition is a condition on uu_t, so we
 * hand apply_bcs_sommerfeld() the interior uu_t = (uu^{n+1} - uu^n)/dt and
 * advance the ghost zones with uu^{n+1} = uu^n + dt*uu_t. This has the same
 * steady state as the RK4 boundary treatment.
 *
 * Before the first step, k_even_gfs must be set equal to y_n_gfs (uu^{-1} = uu^0,
 * i.e., vv = 0 initially, as in the RK4 initial guess).
 */
void conformally_flat_BBH_leapfrog_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt) {
  const paramstruct *restrict params = &griddata->params;
#include "./conformally_flat_BBH_set_Cparameters.h"

  const REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
  REAL *restrict u      = griddata->gridfuncs.y_n_gfs;     // uu^n
  REAL *restrict u_prev = griddata->gridfuncs.k_even_gfs;  // uu^{n-1}, overwritten with uu^{n+1}
  REAL *restrict rhs    = griddata->gridfuncs.k_odd_gfs;
  REAL *restrict residual = griddata->gridfuncs.diagnostic_output_gfs;

  const REAL a = 0.5*eta_damping*dt;
  const REAL inv_1pa = 1.0/(1.0 + a);
  const REAL dt2 = dt*dt;
  const REAL invdt = 1.0/dt;

  // Step 1: rhs[VVGF] = c^2 R(uu^n) on the interior.
  conformally_flat_BBH_rhs_eval(params, &griddata->rfmstruct, auxevol_gfs, u, rhs);

  // Step 2: Interior update. uu^{n+1} goes into u_prev; uu^n stays intact because
  //         the boundary condition below still needs it.
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const int idx = IDX3S(i0,i1,i2);
        const REAL c = auxevol_gfs[IDX4ptS(WAVESPEEDGF,idx)];
        const REAL c2R = rhs[IDX4ptS(VVGF,idx)];
        const REAL un = u[IDX4ptS(UUGF,idx)];
        const REAL unp1 = (2.0*un - (1.0 - a)*u_prev[IDX4ptS(UUGF,idx)] + dt2*c2R)*inv_1pa;
        residual[IDX4ptS(UUGF,idx)] = c2R/(c*c);
        rhs[IDX4ptS(UUGF,idx)] = (unp1 - un)*invdt;
        u_prev[IDX4ptS(UUGF,idx)] = unp1;
      }
    }
  }

  // Step 3: Fill uu_t in the outer (Sommerfeld) and inner (parity) ghost zones.
  conformally_flat_BBH_apply_bcs_sommerfeld(params, griddata->xx, &griddata->bcstruct, 1, evol_gf_parity, u, rhs);

  // Step 4: Advance the ghost zones.
#pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        const int i0i1i2[3] = { i0,i1,i2 };
        if(IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) continue;
        const int idx = IDX3S(i0,i1,i2);
        u_prev[IDX4ptS(UUGF,idx)] = u[IDX4ptS(UUGF,idx)] + dt*rhs[IDX4ptS(UUGF,idx)];
      }
    }
  }

  // Step 5: Rotate time levels: y_n_gfs <- uu^{n+1}, k_even_gfs <- uu^n.
  griddata->gridfuncs.y_n_gfs    = u_prev;
  griddata->gridfuncs.k_even_gfs = u;
}
