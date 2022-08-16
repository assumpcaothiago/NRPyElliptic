#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"
/*
 * Evaluate the RHSs of the Hamiltonian constraint equation
 */
void conformally_flat_BBH_rhs_eval(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "./conformally_flat_BBH_set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    #include "rfm_files/conformally_flat_BBH_rfm_struct__read2.h"
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      #include "rfm_files/conformally_flat_BBH_rfm_struct__read1.h"
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        #include "rfm_files/conformally_flat_BBH_rfm_struct__read0.h"
        {
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          /*
           *  Original SymPy expressions:
           *  "[const double uu_dD0 = invdx0*(-5*uu_i0m1_i1_i2/6 + 5*uu_i0m2_i1_i2/21 - 5*uu_i0m3_i1_i2/84 + 5*uu_i0m4_i1_i2/504 - uu_i0m5_i1_i2/1260 + 5*uu_i0p1_i1_i2/6 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/84 - 5*uu_i0p4_i1_i2/504 + uu_i0p5_i1_i2/1260),
           *    const double uu_dD1 = invdx1*(-5*uu_i0_i1m1_i2/6 + 5*uu_i0_i1m2_i2/21 - 5*uu_i0_i1m3_i2/84 + 5*uu_i0_i1m4_i2/504 - uu_i0_i1m5_i2/1260 + 5*uu_i0_i1p1_i2/6 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/84 - 5*uu_i0_i1p4_i2/504 + uu_i0_i1p5_i2/1260),
           *    const double uu_dDD00 = invdx0**2*(-5269*uu/1800 + 5*uu_i0m1_i1_i2/3 - 5*uu_i0m2_i1_i2/21 + 5*uu_i0m3_i1_i2/126 - 5*uu_i0m4_i1_i2/1008 + uu_i0m5_i1_i2/3150 + 5*uu_i0p1_i1_i2/3 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/126 - 5*uu_i0p4_i1_i2/1008 + uu_i0p5_i1_i2/3150),
           *    const double uu_dDD11 = invdx1**2*(-5269*uu/1800 + 5*uu_i0_i1m1_i2/3 - 5*uu_i0_i1m2_i2/21 + 5*uu_i0_i1m3_i2/126 - 5*uu_i0_i1m4_i2/1008 + uu_i0_i1m5_i2/3150 + 5*uu_i0_i1p1_i2/3 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/126 - 5*uu_i0_i1p4_i2/1008 + uu_i0_i1p5_i2/3150),
           *    const double uu_dDD22 = invdx2**2*(-5269*uu/1800 + 5*uu_i0_i1_i2m1/3 - 5*uu_i0_i1_i2m2/21 + 5*uu_i0_i1_i2m3/126 - 5*uu_i0_i1_i2m4/1008 + uu_i0_i1_i2m5/3150 + 5*uu_i0_i1_i2p1/3 - 5*uu_i0_i1_i2p2/21 + 5*uu_i0_i1_i2p3/126 - 5*uu_i0_i1_i2p4/1008 + uu_i0_i1_i2p5/3150)]"
           */
          const double wavespeed = auxevol_gfs[IDX4S(WAVESPEEDGF, i0,i1,i2)];
          const double psi_background = auxevol_gfs[IDX4S(PSI_BACKGROUNDGF, i0,i1,i2)];
          const double ADD_times_AUU = auxevol_gfs[IDX4S(ADD_TIMES_AUUGF, i0,i1,i2)];
          const double uu_i0_i1_i2m5 = in_gfs[IDX4S(UUGF, i0,i1,i2-5)];
          const double uu_i0_i1_i2m4 = in_gfs[IDX4S(UUGF, i0,i1,i2-4)];
          const double uu_i0_i1_i2m3 = in_gfs[IDX4S(UUGF, i0,i1,i2-3)];
          const double uu_i0_i1_i2m2 = in_gfs[IDX4S(UUGF, i0,i1,i2-2)];
          const double uu_i0_i1_i2m1 = in_gfs[IDX4S(UUGF, i0,i1,i2-1)];
          const double uu_i0_i1m5_i2 = in_gfs[IDX4S(UUGF, i0,i1-5,i2)];
          const double uu_i0_i1m4_i2 = in_gfs[IDX4S(UUGF, i0,i1-4,i2)];
          const double uu_i0_i1m3_i2 = in_gfs[IDX4S(UUGF, i0,i1-3,i2)];
          const double uu_i0_i1m2_i2 = in_gfs[IDX4S(UUGF, i0,i1-2,i2)];
          const double uu_i0_i1m1_i2 = in_gfs[IDX4S(UUGF, i0,i1-1,i2)];
          const double uu_i0m5_i1_i2 = in_gfs[IDX4S(UUGF, i0-5,i1,i2)];
          const double uu_i0m4_i1_i2 = in_gfs[IDX4S(UUGF, i0-4,i1,i2)];
          const double uu_i0m3_i1_i2 = in_gfs[IDX4S(UUGF, i0-3,i1,i2)];
          const double uu_i0m2_i1_i2 = in_gfs[IDX4S(UUGF, i0-2,i1,i2)];
          const double uu_i0m1_i1_i2 = in_gfs[IDX4S(UUGF, i0-1,i1,i2)];
          const double uu = in_gfs[IDX4S(UUGF, i0,i1,i2)];
          const double uu_i0p1_i1_i2 = in_gfs[IDX4S(UUGF, i0+1,i1,i2)];
          const double uu_i0p2_i1_i2 = in_gfs[IDX4S(UUGF, i0+2,i1,i2)];
          const double uu_i0p3_i1_i2 = in_gfs[IDX4S(UUGF, i0+3,i1,i2)];
          const double uu_i0p4_i1_i2 = in_gfs[IDX4S(UUGF, i0+4,i1,i2)];
          const double uu_i0p5_i1_i2 = in_gfs[IDX4S(UUGF, i0+5,i1,i2)];
          const double uu_i0_i1p1_i2 = in_gfs[IDX4S(UUGF, i0,i1+1,i2)];
          const double uu_i0_i1p2_i2 = in_gfs[IDX4S(UUGF, i0,i1+2,i2)];
          const double uu_i0_i1p3_i2 = in_gfs[IDX4S(UUGF, i0,i1+3,i2)];
          const double uu_i0_i1p4_i2 = in_gfs[IDX4S(UUGF, i0,i1+4,i2)];
          const double uu_i0_i1p5_i2 = in_gfs[IDX4S(UUGF, i0,i1+5,i2)];
          const double uu_i0_i1_i2p1 = in_gfs[IDX4S(UUGF, i0,i1,i2+1)];
          const double uu_i0_i1_i2p2 = in_gfs[IDX4S(UUGF, i0,i1,i2+2)];
          const double uu_i0_i1_i2p3 = in_gfs[IDX4S(UUGF, i0,i1,i2+3)];
          const double uu_i0_i1_i2p4 = in_gfs[IDX4S(UUGF, i0,i1,i2+4)];
          const double uu_i0_i1_i2p5 = in_gfs[IDX4S(UUGF, i0,i1,i2+5)];
          const double vv = in_gfs[IDX4S(VVGF, i0,i1,i2)];
          const double FDPart1_Rational_5_6 = 5.0/6.0;
          const double FDPart1_Rational_5_21 = 5.0/21.0;
          const double FDPart1_Rational_5_84 = 5.0/84.0;
          const double FDPart1_Rational_5_504 = 5.0/504.0;
          const double FDPart1_Rational_1_1260 = 1.0/1260.0;
          const double FDPart1_Rational_5269_1800 = 5269.0/1800.0;
          const double FDPart1_Rational_5_1008 = 5.0/1008.0;
          const double FDPart1_Rational_1_3150 = 1.0/3150.0;
          const double FDPart1_Rational_5_3 = 5.0/3.0;
          const double FDPart1_Rational_5_126 = 5.0/126.0;
          const double FDPart1_0 = -FDPart1_Rational_5269_1800*uu;
          const double uu_dD0 = invdx0*(FDPart1_Rational_1_1260*(-uu_i0m5_i1_i2 + uu_i0p5_i1_i2) + FDPart1_Rational_5_21*(uu_i0m2_i1_i2 - uu_i0p2_i1_i2) + FDPart1_Rational_5_504*(uu_i0m4_i1_i2 - uu_i0p4_i1_i2) + FDPart1_Rational_5_6*(-uu_i0m1_i1_i2 + uu_i0p1_i1_i2) + FDPart1_Rational_5_84*(-uu_i0m3_i1_i2 + uu_i0p3_i1_i2));
          const double uu_dD1 = invdx1*(FDPart1_Rational_1_1260*(-uu_i0_i1m5_i2 + uu_i0_i1p5_i2) + FDPart1_Rational_5_21*(uu_i0_i1m2_i2 - uu_i0_i1p2_i2) + FDPart1_Rational_5_504*(uu_i0_i1m4_i2 - uu_i0_i1p4_i2) + FDPart1_Rational_5_6*(-uu_i0_i1m1_i2 + uu_i0_i1p1_i2) + FDPart1_Rational_5_84*(-uu_i0_i1m3_i2 + uu_i0_i1p3_i2));
          const double uu_dDD00 = ((invdx0)*(invdx0))*(FDPart1_0 + FDPart1_Rational_1_3150*(uu_i0m5_i1_i2 + uu_i0p5_i1_i2) + FDPart1_Rational_5_1008*(-uu_i0m4_i1_i2 - uu_i0p4_i1_i2) + FDPart1_Rational_5_126*(uu_i0m3_i1_i2 + uu_i0p3_i1_i2) + FDPart1_Rational_5_21*(-uu_i0m2_i1_i2 - uu_i0p2_i1_i2) + FDPart1_Rational_5_3*(uu_i0m1_i1_i2 + uu_i0p1_i1_i2));
          const double uu_dDD11 = ((invdx1)*(invdx1))*(FDPart1_0 + FDPart1_Rational_1_3150*(uu_i0_i1m5_i2 + uu_i0_i1p5_i2) + FDPart1_Rational_5_1008*(-uu_i0_i1m4_i2 - uu_i0_i1p4_i2) + FDPart1_Rational_5_126*(uu_i0_i1m3_i2 + uu_i0_i1p3_i2) + FDPart1_Rational_5_21*(-uu_i0_i1m2_i2 - uu_i0_i1p2_i2) + FDPart1_Rational_5_3*(uu_i0_i1m1_i2 + uu_i0_i1p1_i2));
          const double uu_dDD22 = ((invdx2)*(invdx2))*(FDPart1_0 + FDPart1_Rational_1_3150*(uu_i0_i1_i2m5 + uu_i0_i1_i2p5) + FDPart1_Rational_5_1008*(-uu_i0_i1_i2m4 - uu_i0_i1_i2p4) + FDPart1_Rational_5_126*(uu_i0_i1_i2m3 + uu_i0_i1_i2p3) + FDPart1_Rational_5_21*(-uu_i0_i1_i2m2 - uu_i0_i1_i2p2) + FDPart1_Rational_5_3*(uu_i0_i1_i2m1 + uu_i0_i1_i2p1));
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          /*
           *  Original SymPy expressions:
           *  "[rhs_gfs[IDX4S(UUGF, i0, i1, i2)] = -eta_damping*uu + vv,
           *    rhs_gfs[IDX4S(VVGF, i0, i1, i2)] = wavespeed**2*(ADD_times_AUU/(8*(psi_background + uu)**7) - uu_dD0*(-f2_of_xx0_xx1__D0*f3_of_xx0**2/(f0_of_xx0__D0**2*f2_of_xx0_xx1**3) + f3_of_xx0*(-f0_of_xx0__D0*f2_of_xx0_xx1*f3_of_xx0__D0 + f3_of_xx0*(f0_of_xx0__D0*f2_of_xx0_xx1__D0 + f0_of_xx0__DD00*f2_of_xx0_xx1))/(f0_of_xx0__D0**3*f2_of_xx0_xx1**3) - f3_of_xx0**2/(f0_of_xx0*f0_of_xx0__D0*f2_of_xx0_xx1**2)) + uu_dDD11/f2_of_xx0_xx1**2 + f1_of_xx1__D1*uu_dD1/(f1_of_xx1*f2_of_xx0_xx1**2) + f3_of_xx0**2*uu_dDD00/(f0_of_xx0__D0**2*f2_of_xx0_xx1**2) + uu_dDD22/(f0_of_xx0**2*f1_of_xx1**2))]"
           */
          const double FDPart3_0 = (1.0/((f2_of_xx0_xx1)*(f2_of_xx0_xx1)));
          const double FDPart3_2 = ((f3_of_xx0)*(f3_of_xx0))/((f0_of_xx0__D0)*(f0_of_xx0__D0));
          const double FDPart3_3 = (1.0/((f2_of_xx0_xx1)*(f2_of_xx0_xx1)*(f2_of_xx0_xx1)));
          rhs_gfs[IDX4S(UUGF, i0, i1, i2)] = -eta_damping*uu + vv;
          rhs_gfs[IDX4S(VVGF, i0, i1, i2)] = ((wavespeed)*(wavespeed))*((1.0/8.0)*ADD_times_AUU/pow(psi_background + uu, 7) + FDPart3_0*FDPart3_2*uu_dDD00 + FDPart3_0*uu_dDD11 + FDPart3_0*f1_of_xx1__D1*uu_dD1/f1_of_xx1 - uu_dD0*(-FDPart3_0*((f3_of_xx0)*(f3_of_xx0))/(f0_of_xx0*f0_of_xx0__D0) - FDPart3_2*FDPart3_3*f2_of_xx0_xx1__D0 + FDPart3_3*f3_of_xx0*(-f0_of_xx0__D0*f2_of_xx0_xx1*f3_of_xx0__D0 + f3_of_xx0*(f0_of_xx0__D0*f2_of_xx0_xx1__D0 + f0_of_xx0__DD00*f2_of_xx0_xx1))/((f0_of_xx0__D0)*(f0_of_xx0__D0)*(f0_of_xx0__D0))) + uu_dDD22/(((f0_of_xx0)*(f0_of_xx0))*((f1_of_xx1)*(f1_of_xx1))));
        }
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
