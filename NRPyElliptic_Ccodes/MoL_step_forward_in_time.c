#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "apply_bcs_sommerfeld.h"
#include "./SIMD/SIMD_intrinsics.h"
/*
 * Method of Lines (MoL) for "RK4" method: Step forward one full timestep.
 */
void MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt) {

  // C code implementation of -={ RK4 }=- Method of Lines timestepping.

  // -={ START k1 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    const rfm_struct *restrict rfmstruct = &griddata->rfmstruct;
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    rhs_eval(params, rfmstruct, auxevol_gfs, y_n_gfs, k_odd_gfs);apply_bcs_sommerfeld(params, griddata->xx, bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs, k_odd_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=SIMD_width) {
      const REAL_SIMD_ARRAY k_odd_gfsL = ReadSIMD(&k_odd_gfs[i]);
      const REAL_SIMD_ARRAY y_n_gfsL = ReadSIMD(&y_n_gfs[i]);
      const REAL_SIMD_ARRAY DT = ConstSIMD(dt);
      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

        const REAL_SIMD_ARRAY __RHS_exp_0 = MulSIMD(_Rational_1_6, MulSIMD(DT, k_odd_gfsL));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(_Rational_1_2, MulSIMD(DT, k_odd_gfsL), y_n_gfsL);
    WriteSIMD(&y_nplus1_running_total_gfs[i], __RHS_exp_0);
    WriteSIMD(&k_odd_gfs[i], __RHS_exp_1);
    }
  }
  // -={ END k1 substep }=-

  // -={ START k2 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    const rfm_struct *restrict rfmstruct = &griddata->rfmstruct;
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    rhs_eval(params, rfmstruct, auxevol_gfs, k_odd_gfs, k_even_gfs);apply_bcs_sommerfeld(params, griddata->xx, bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs, k_even_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=SIMD_width) {
      const REAL_SIMD_ARRAY k_even_gfsL = ReadSIMD(&k_even_gfs[i]);
      const REAL_SIMD_ARRAY y_nplus1_running_total_gfsL = ReadSIMD(&y_nplus1_running_total_gfs[i]);
      const REAL_SIMD_ARRAY y_n_gfsL = ReadSIMD(&y_n_gfs[i]);
      const REAL_SIMD_ARRAY DT = ConstSIMD(dt);
      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_3 = 1.0/3.0;
      const REAL_SIMD_ARRAY _Rational_1_3 = ConstSIMD(tmp_Rational_1_3);

        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(_Rational_1_3, MulSIMD(DT, k_even_gfsL), y_nplus1_running_total_gfsL);
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(_Rational_1_2, MulSIMD(DT, k_even_gfsL), y_n_gfsL);
    WriteSIMD(&y_nplus1_running_total_gfs[i], __RHS_exp_0);
    WriteSIMD(&k_even_gfs[i], __RHS_exp_1);
    }
  }
  // -={ END k2 substep }=-

  // -={ START k3 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    const rfm_struct *restrict rfmstruct = &griddata->rfmstruct;
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    rhs_eval(params, rfmstruct, auxevol_gfs, k_even_gfs, k_odd_gfs);apply_bcs_sommerfeld(params, griddata->xx, bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_even_gfs, k_odd_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=SIMD_width) {
      const REAL_SIMD_ARRAY k_odd_gfsL = ReadSIMD(&k_odd_gfs[i]);
      const REAL_SIMD_ARRAY y_nplus1_running_total_gfsL = ReadSIMD(&y_nplus1_running_total_gfs[i]);
      const REAL_SIMD_ARRAY y_n_gfsL = ReadSIMD(&y_n_gfs[i]);
      const REAL_SIMD_ARRAY DT = ConstSIMD(dt);
      const double tmp_Rational_1_3 = 1.0/3.0;
      const REAL_SIMD_ARRAY _Rational_1_3 = ConstSIMD(tmp_Rational_1_3);

        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(_Rational_1_3, MulSIMD(DT, k_odd_gfsL), y_nplus1_running_total_gfsL);
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(DT, k_odd_gfsL, y_n_gfsL);
    WriteSIMD(&y_nplus1_running_total_gfs[i], __RHS_exp_0);
    WriteSIMD(&k_odd_gfs[i], __RHS_exp_1);
    }
  }
  // -={ END k3 substep }=-

  // -={ START k4 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    const rfm_struct *restrict rfmstruct = &griddata->rfmstruct;
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    rhs_eval(params, rfmstruct, auxevol_gfs, k_odd_gfs, k_even_gfs);apply_bcs_sommerfeld(params, griddata->xx, bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs, k_even_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=SIMD_width) {
      const REAL_SIMD_ARRAY k_even_gfsL = ReadSIMD(&k_even_gfs[i]);
      const REAL_SIMD_ARRAY y_n_gfsL = ReadSIMD(&y_n_gfs[i]);
      const REAL_SIMD_ARRAY y_nplus1_running_total_gfsL = ReadSIMD(&y_nplus1_running_total_gfs[i]);
      const REAL_SIMD_ARRAY DT = ConstSIMD(dt);
      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

        const REAL_SIMD_ARRAY __RHS_exp_0 = AddSIMD(y_n_gfsL, FusedMulAddSIMD(_Rational_1_6, MulSIMD(DT, k_even_gfsL), y_nplus1_running_total_gfsL));
    WriteSIMD(&y_n_gfs[i], __RHS_exp_0);
    }
  }
  // -={ END k4 substep }=-

}
