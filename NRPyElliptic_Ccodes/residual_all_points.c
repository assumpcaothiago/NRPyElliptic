#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "./SIMD/SIMD_intrinsics.h"
/*
 * Evaluate the residual at all points
 */
void residual_all_points(const paramstruct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict residual_gf) {
#include "./set_Cparameters-SIMD.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    #include "rfm_files/rfm_struct__SIMD_outer_read2.h"
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      #include "rfm_files/rfm_struct__SIMD_outer_read1.h"
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += SIMD_width) {
        #include "rfm_files/rfm_struct__SIMD_inner_read0.h"
        {
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          /*
           *  Original SymPy expressions:
           *  "[const REAL_SIMD_ARRAY uu_dD0 = invdx0*(-5*uu_i0m1_i1_i2/6 + 5*uu_i0m2_i1_i2/21 - 5*uu_i0m3_i1_i2/84 + 5*uu_i0m4_i1_i2/504 - uu_i0m5_i1_i2/1260 + 5*uu_i0p1_i1_i2/6 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/84 - 5*uu_i0p4_i1_i2/504 + uu_i0p5_i1_i2/1260),
           *    const REAL_SIMD_ARRAY uu_dD1 = invdx1*(-5*uu_i0_i1m1_i2/6 + 5*uu_i0_i1m2_i2/21 - 5*uu_i0_i1m3_i2/84 + 5*uu_i0_i1m4_i2/504 - uu_i0_i1m5_i2/1260 + 5*uu_i0_i1p1_i2/6 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/84 - 5*uu_i0_i1p4_i2/504 + uu_i0_i1p5_i2/1260),
           *    const REAL_SIMD_ARRAY uu_dDD00 = invdx0**2*(-5269*uu/1800 + 5*uu_i0m1_i1_i2/3 - 5*uu_i0m2_i1_i2/21 + 5*uu_i0m3_i1_i2/126 - 5*uu_i0m4_i1_i2/1008 + uu_i0m5_i1_i2/3150 + 5*uu_i0p1_i1_i2/3 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/126 - 5*uu_i0p4_i1_i2/1008 + uu_i0p5_i1_i2/3150),
           *    const REAL_SIMD_ARRAY uu_dDD11 = invdx1**2*(-5269*uu/1800 + 5*uu_i0_i1m1_i2/3 - 5*uu_i0_i1m2_i2/21 + 5*uu_i0_i1m3_i2/126 - 5*uu_i0_i1m4_i2/1008 + uu_i0_i1m5_i2/3150 + 5*uu_i0_i1p1_i2/3 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/126 - 5*uu_i0_i1p4_i2/1008 + uu_i0_i1p5_i2/3150),
           *    const REAL_SIMD_ARRAY uu_dDD22 = invdx2**2*(-5269*uu/1800 + 5*uu_i0_i1_i2m1/3 - 5*uu_i0_i1_i2m2/21 + 5*uu_i0_i1_i2m3/126 - 5*uu_i0_i1_i2m4/1008 + uu_i0_i1_i2m5/3150 + 5*uu_i0_i1_i2p1/3 - 5*uu_i0_i1_i2p2/21 + 5*uu_i0_i1_i2p3/126 - 5*uu_i0_i1_i2p4/1008 + uu_i0_i1_i2p5/3150)]"
           */
          const REAL_SIMD_ARRAY psi_background = ReadSIMD(&auxevol_gfs[IDX4S(PSI_BACKGROUNDGF, i0,i1,i2)]);
          const REAL_SIMD_ARRAY ADD_times_AUU = ReadSIMD(&auxevol_gfs[IDX4S(ADD_TIMES_AUUGF, i0,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2m5 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-5)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2m4 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-4)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2m3 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-3)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-2)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-1)]);
          const REAL_SIMD_ARRAY uu_i0_i1m5_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-5,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1m4_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-4,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1m3_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-3,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-2,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-1,i2)]);
          const REAL_SIMD_ARRAY uu_i0m5_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-5,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0m4_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-4,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0m3_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-3,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-2,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-1,i1,i2)]);
          const REAL_SIMD_ARRAY uu = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+1,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+2,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0p3_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+3,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0p4_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+4,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0p5_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+5,i1,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+1,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+2,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1p3_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+3,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1p4_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+4,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1p5_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+5,i2)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+1)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+2)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2p3 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+3)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2p4 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+4)]);
          const REAL_SIMD_ARRAY uu_i0_i1_i2p5 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+5)]);
          const double tmpFDPart1_NegativeOne_ = -1.0;
          const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(tmpFDPart1_NegativeOne_);
          const double tmpFDPart1_Rational_1_1260 = 1.0/1260.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_1_1260 = ConstSIMD(tmpFDPart1_Rational_1_1260);
          const double tmpFDPart1_Rational_1_3150 = 1.0/3150.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_1_3150 = ConstSIMD(tmpFDPart1_Rational_1_3150);
          const double tmpFDPart1_Rational_5269_1800 = 5269.0/1800.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5269_1800 = ConstSIMD(tmpFDPart1_Rational_5269_1800);
          const double tmpFDPart1_Rational_5_1008 = 5.0/1008.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_1008 = ConstSIMD(tmpFDPart1_Rational_5_1008);
          const double tmpFDPart1_Rational_5_126 = 5.0/126.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_126 = ConstSIMD(tmpFDPart1_Rational_5_126);
          const double tmpFDPart1_Rational_5_21 = 5.0/21.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_21 = ConstSIMD(tmpFDPart1_Rational_5_21);
          const double tmpFDPart1_Rational_5_3 = 5.0/3.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_3 = ConstSIMD(tmpFDPart1_Rational_5_3);
          const double tmpFDPart1_Rational_5_504 = 5.0/504.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_504 = ConstSIMD(tmpFDPart1_Rational_5_504);
          const double tmpFDPart1_Rational_5_6 = 5.0/6.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(tmpFDPart1_Rational_5_6);
          const double tmpFDPart1_Rational_5_84 = 5.0/84.0;
          const REAL_SIMD_ARRAY FDPart1_Rational_5_84 = ConstSIMD(tmpFDPart1_Rational_5_84);
          const REAL_SIMD_ARRAY FDPart1_0 = MulSIMD(FDPart1_Rational_5269_1800, uu);
          const REAL_SIMD_ARRAY uu_dD0 = MulSIMD(invdx0, FusedMulAddSIMD(FDPart1_Rational_5_504, SubSIMD(uu_i0m4_i1_i2, uu_i0p4_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_6, SubSIMD(uu_i0p1_i1_i2, uu_i0m1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_84, SubSIMD(uu_i0p3_i1_i2, uu_i0m3_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_1260, SubSIMD(uu_i0p5_i1_i2, uu_i0m5_i1_i2), MulSIMD(FDPart1_Rational_5_21, SubSIMD(uu_i0m2_i1_i2, uu_i0p2_i1_i2)))))));
          const REAL_SIMD_ARRAY uu_dD1 = MulSIMD(invdx1, FusedMulAddSIMD(FDPart1_Rational_5_504, SubSIMD(uu_i0_i1m4_i2, uu_i0_i1p4_i2), FusedMulAddSIMD(FDPart1_Rational_5_6, SubSIMD(uu_i0_i1p1_i2, uu_i0_i1m1_i2), FusedMulAddSIMD(FDPart1_Rational_5_84, SubSIMD(uu_i0_i1p3_i2, uu_i0_i1m3_i2), FusedMulAddSIMD(FDPart1_Rational_1_1260, SubSIMD(uu_i0_i1p5_i2, uu_i0_i1m5_i2), MulSIMD(FDPart1_Rational_5_21, SubSIMD(uu_i0_i1m2_i2, uu_i0_i1p2_i2)))))));
          const REAL_SIMD_ARRAY uu_dDD00 = MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0m3_i1_i2, uu_i0p3_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0m1_i1_i2, uu_i0p1_i1_i2), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0m5_i1_i2, uu_i0p5_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0m4_i1_i2, uu_i0p4_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0m2_i1_i2, uu_i0p2_i1_i2), FDPart1_0))))));
          const REAL_SIMD_ARRAY uu_dDD11 = MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0_i1m3_i2, uu_i0_i1p3_i2), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0_i1m1_i2, uu_i0_i1p1_i2), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0_i1m5_i2, uu_i0_i1p5_i2), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0_i1m4_i2, uu_i0_i1p4_i2), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0_i1m2_i2, uu_i0_i1p2_i2), FDPart1_0))))));
          const REAL_SIMD_ARRAY uu_dDD22 = MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0_i1_i2m3, uu_i0_i1_i2p3), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0_i1_i2m1, uu_i0_i1_i2p1), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0_i1_i2m5, uu_i0_i1_i2p5), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0_i1_i2m4, uu_i0_i1_i2p4), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0_i1_i2m2, uu_i0_i1_i2p2), FDPart1_0))))));
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          /*
           *  Original SymPy expression:
           *  "const REAL_SIMD_ARRAY __RHS_exp_0 = ADD_times_AUU/(8*(psi_background + uu)**7) - uu_dD0*(-f2_of_xx0_xx1__D0*f3_of_xx0**2/(f0_of_xx0__D0**2*f2_of_xx0_xx1**3) + f3_of_xx0*(-f0_of_xx0__D0*f2_of_xx0_xx1*f3_of_xx0__D0 + f3_of_xx0*(f0_of_xx0__D0*f2_of_xx0_xx1__D0 + f0_of_xx0__DD00*f2_of_xx0_xx1))/(f0_of_xx0__D0**3*f2_of_xx0_xx1**3) - f3_of_xx0**2/(f0_of_xx0*f0_of_xx0__D0*f2_of_xx0_xx1**2)) + uu_dDD11/f2_of_xx0_xx1**2 + f1_of_xx1__D1*uu_dD1/(f1_of_xx1*f2_of_xx0_xx1**2) + f3_of_xx0**2*uu_dDD00/(f0_of_xx0__D0**2*f2_of_xx0_xx1**2) + uu_dDD22/(f0_of_xx0**2*f1_of_xx1**2)"
           */
          const double tmpFDPart3_Integer_1 = 1.0;
          const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(tmpFDPart3_Integer_1);
          const double tmpFDPart3_NegativeOne_ = -1.0;
          const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(tmpFDPart3_NegativeOne_);
          const double tmpFDPart3_Rational_1_8 = 1.0/8.0;
          const REAL_SIMD_ARRAY FDPart3_Rational_1_8 = ConstSIMD(tmpFDPart3_Rational_1_8);
          const REAL_SIMD_ARRAY FDPart3_0 = DivSIMD(FDPart3_Integer_1, MulSIMD(f2_of_xx0_xx1, f2_of_xx0_xx1));
          const REAL_SIMD_ARRAY FDPart3_2 = DivSIMD(MulSIMD(f3_of_xx0, f3_of_xx0), MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0));
          const REAL_SIMD_ARRAY FDPart3_3 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(f2_of_xx0_xx1, f2_of_xx0_xx1), f2_of_xx0_xx1));
          const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(uu_dDD22, DivSIMD(DivSIMD(FDPart3_Integer_1, MulSIMD(f1_of_xx1, f1_of_xx1)), MulSIMD(f0_of_xx0, f0_of_xx0)), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_2, uu_dDD00), FusedMulAddSIMD(FDPart3_Rational_1_8, DivSIMD(ADD_times_AUU, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(AddSIMD(psi_background, uu), AddSIMD(psi_background, uu)), AddSIMD(psi_background, uu)), AddSIMD(psi_background, uu)), AddSIMD(psi_background, uu)), AddSIMD(psi_background, uu)), AddSIMD(psi_background, uu))), FusedMulAddSIMD(MulSIMD(FDPart3_0, f1_of_xx1__D1), DivSIMD(uu_dD1, f1_of_xx1), FusedMulSubSIMD(FDPart3_0, uu_dDD11, MulSIMD(uu_dD0, FusedMulSubSIMD(MulSIMD(FDPart3_3, f3_of_xx0), DivSIMD(FusedMulSubSIMD(f3_of_xx0, FusedMulAddSIMD(f0_of_xx0__D0, f2_of_xx0_xx1__D0, MulSIMD(f0_of_xx0__DD00, f2_of_xx0_xx1)), MulSIMD(f0_of_xx0__D0, MulSIMD(f2_of_xx0_xx1, f3_of_xx0__D0))), MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0)), FusedMulSubSIMD(FDPart3_2, MulSIMD(FDPart3_3, f2_of_xx0_xx1__D0), DivSIMD(MulSIMD(MulSIMD(FDPart3_0, FDPart3_NegativeOne_), DivSIMD(MulSIMD(f3_of_xx0, f3_of_xx0), f0_of_xx0__D0)), f0_of_xx0)))))))));
          WriteSIMD(&residual_gf[IDX4S(UUGF, i0, i1, i2)], __RHS_exp_0);
        }
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += SIMD_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
