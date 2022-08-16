#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"

/*
 * This function sets ADM quantities from the initial data variables
 */
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
                                                                   REAL *restrict KDD22) {
#include "./conformally_flat_BBH_set_Cparameters.h"

    if( CCTK_EQUALS(initial_lapse,"NRPyEllipticET-psi^n") ) {
      conformally_flat_BBH_lapse_power_of_psi_initial_data(cctkGH,params,x,y,z,in_gfs,alpha);
    }
    else if( CCTK_EQUALS(initial_lapse,"NRPyEllipticET-averaged") ) {
      conformally_flat_BBH_lapse_averaged_initial_data(    cctkGH,orbital_plane,params,x,y,z,position_shift,in_gfs,alpha);
    }
    else {
      CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                  "NRPyEllipticET does not support the lapse initial data: %s",initial_lapse);
    }

    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {

                const int idx      = CCTK_GFINDEX3D(cctkGH,i0,i1,i2);

                const double Cartx = x[idx];
                const double Carty = y[idx];
                const double Cartz = z[idx];
                /*
                 * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                 */
                const double uu = in_gfs[IDX4S(UUGF, i0,i1,i2)];
                /*
                 * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                 */
                const double FDPart3_1 = puncture_0_x - Cartx;
                const double FDPart3_4 = puncture_0_y - Carty;
                const double FDPart3_7 = puncture_0_z - Cartz;
                const double FDPart3_9 = ((FDPart3_1)*(FDPart3_1)) + ((FDPart3_4)*(FDPart3_4)) + ((FDPart3_7)*(FDPart3_7));
                const double FDPart3_10 = puncture_1_x - Cartx;
                const double FDPart3_12 = puncture_1_y - Carty;
                const double FDPart3_14 = puncture_1_z - Cartz;
                const double FDPart3_16 = ((FDPart3_10)*(FDPart3_10)) + ((FDPart3_12)*(FDPart3_12)) + ((FDPart3_14)*(FDPart3_14));
                const double FDPart3_17 = uu + 1 + (1.0/2.0)*puncture_0_bare_mass/sqrt(FDPart3_9) + (1.0/2.0)*puncture_1_bare_mass/sqrt(FDPart3_16);
                const double FDPart3_18 = ((FDPart3_17)*(FDPart3_17)*(FDPart3_17)*(FDPart3_17));
                const double FDPart3_19 = (1.0/((FDPart3_17)*(FDPart3_17)));
                const double FDPart3_20 = -FDPart3_1*puncture_0_P_x;
                const double FDPart3_21 = pow(FDPart3_9, -3.0/2.0);
                const double FDPart3_22 = (1.0/3.0)*FDPart3_21;
                const double FDPart3_23 = FDPart3_21*puncture_0_P_x;
                const double FDPart3_25 = -FDPart3_10*puncture_1_P_x;
                const double FDPart3_26 = pow(FDPart3_16, -3.0/2.0);
                const double FDPart3_27 = (1.0/3.0)*FDPart3_26;
                const double FDPart3_28 = FDPart3_26*puncture_1_P_x;
                const double FDPart3_31 = 3*puncture_0_x - 3*Cartx;
                const double FDPart3_32 = -FDPart3_4*puncture_0_S_z + FDPart3_7*puncture_0_S_y;
                const double FDPart3_33 = pow(FDPart3_9, -5.0/2.0);
                const double FDPart3_34 = FDPart3_1*puncture_0_S_z - FDPart3_7*puncture_0_S_x;
                const double FDPart3_36 = 3*puncture_0_y - 3*Carty;
                const double FDPart3_38 = FDPart3_33*FDPart3_34*FDPart3_36;
                const double FDPart3_40 = -FDPart3_1*puncture_0_S_y + FDPart3_4*puncture_0_S_x;
                const double FDPart3_42 = 3*puncture_0_z - 3*Cartz;
                const double FDPart3_44 = FDPart3_33*FDPart3_40*FDPart3_42;
                const double FDPart3_46 = 3*puncture_1_x - 3*Cartx;
                const double FDPart3_47 = -FDPart3_12*puncture_1_S_z + FDPart3_14*puncture_1_S_y;
                const double FDPart3_48 = pow(FDPart3_16, -5.0/2.0);
                const double FDPart3_49 = FDPart3_10*puncture_1_S_z - FDPart3_14*puncture_1_S_x;
                const double FDPart3_50 = 3*puncture_1_y - 3*Carty;
                const double FDPart3_52 = FDPart3_48*FDPart3_49*FDPart3_50;
                const double FDPart3_54 = -FDPart3_10*puncture_1_S_y + FDPart3_12*puncture_1_S_x;
                const double FDPart3_55 = 3*puncture_1_z - 3*Cartz;
                const double FDPart3_57 = FDPart3_48*FDPart3_54*FDPart3_55;
                const double FDPart3_59 = -FDPart3_4*puncture_0_P_y;
                const double FDPart3_60 = -FDPart3_7*puncture_0_P_z;
                const double FDPart3_61 = -1.0/4.0*FDPart3_20 - 1.0/4.0*FDPart3_59 - 1.0/4.0*FDPart3_60;
                const double FDPart3_62 = -FDPart3_33*FDPart3_36*FDPart3_4*FDPart3_61;
                const double FDPart3_64 = -FDPart3_33*FDPart3_42*FDPart3_61*FDPart3_7;
                const double FDPart3_66 = -FDPart3_12*puncture_1_P_y;
                const double FDPart3_67 = -FDPart3_14*puncture_1_P_z;
                const double FDPart3_68 = -1.0/4.0*FDPart3_25 - 1.0/4.0*FDPart3_66 - 1.0/4.0*FDPart3_67;
                const double FDPart3_69 = -FDPart3_12*FDPart3_48*FDPart3_50*FDPart3_68;
                const double FDPart3_71 = -FDPart3_14*FDPart3_48*FDPart3_55*FDPart3_68;
                const double FDPart3_75 = FDPart3_21*puncture_0_P_y;
                const double FDPart3_80 = FDPart3_26*puncture_1_P_y;
                const double FDPart3_89 = (1.0/6.0)*FDPart3_21;
                const double FDPart3_90 = (1.0/6.0)*FDPart3_26;
                const double FDPart3_97 = (7.0/6.0)*FDPart3_1*FDPart3_23 + (2.0/3.0)*FDPart3_1*FDPart3_31*FDPart3_33*FDPart3_61 + (7.0/6.0)*FDPart3_10*FDPart3_28 + (2.0/3.0)*FDPart3_10*FDPart3_46*FDPart3_48*FDPart3_68 + FDPart3_20*FDPart3_89 + FDPart3_25*FDPart3_90 - 2.0/3.0*FDPart3_31*FDPart3_32*FDPart3_33 - 2.0/3.0*FDPart3_46*FDPart3_47*FDPart3_48;
                betaU0[idx] = 0;
                betaU1[idx] = 0;
                betaU2[idx] = 0;
                gammaDD00[idx] = FDPart3_18;
                gammaDD01[idx] = 0;
                gammaDD02[idx] = 0;
                gammaDD11[idx] = FDPart3_18;
                gammaDD12[idx] = 0;
                gammaDD22[idx] = FDPart3_18;
                KDD00[idx] = FDPart3_19*(-7.0/3.0*FDPart3_1*FDPart3_23 - 4.0/3.0*FDPart3_1*FDPart3_31*FDPart3_33*FDPart3_61 - 7.0/3.0*FDPart3_10*FDPart3_28 - 4.0/3.0*FDPart3_10*FDPart3_46*FDPart3_48*FDPart3_68 + FDPart3_12*FDPart3_26*puncture_1_P_y + FDPart3_14*FDPart3_26*puncture_1_P_z - FDPart3_20*FDPart3_22 + FDPart3_21*FDPart3_4*puncture_0_P_y + FDPart3_21*FDPart3_7*puncture_0_P_z - FDPart3_25*FDPart3_27 + (4.0/3.0)*FDPart3_31*FDPart3_32*FDPart3_33 - 2.0/3.0*FDPart3_38 - 2.0/3.0*FDPart3_44 + (4.0/3.0)*FDPart3_46*FDPart3_47*FDPart3_48 - 2.0/3.0*FDPart3_52 - 2.0/3.0*FDPart3_57 - 2.0/3.0*FDPart3_62 - 2.0/3.0*FDPart3_64 - 2.0/3.0*FDPart3_69 - 2.0/3.0*FDPart3_71);
                KDD01[idx] = FDPart3_19*(-FDPart3_1*FDPart3_33*FDPart3_36*FDPart3_61 - 3.0/2.0*FDPart3_1*FDPart3_75 - FDPart3_10*FDPart3_48*FDPart3_50*FDPart3_68 - 3.0/2.0*FDPart3_10*FDPart3_80 - 3.0/2.0*FDPart3_12*FDPart3_28 - FDPart3_12*FDPart3_46*FDPart3_48*FDPart3_68 - 3.0/2.0*FDPart3_23*FDPart3_4 + FDPart3_31*FDPart3_33*FDPart3_34 - FDPart3_31*FDPart3_33*FDPart3_4*FDPart3_61 + FDPart3_32*FDPart3_33*FDPart3_36 + FDPart3_46*FDPart3_48*FDPart3_49 + FDPart3_47*FDPart3_48*FDPart3_50);
                KDD02[idx] = FDPart3_19*(-3.0/2.0*FDPart3_1*FDPart3_21*puncture_0_P_z - FDPart3_1*FDPart3_33*FDPart3_42*FDPart3_61 - 3.0/2.0*FDPart3_10*FDPart3_26*puncture_1_P_z - FDPart3_10*FDPart3_48*FDPart3_55*FDPart3_68 - 3.0/2.0*FDPart3_14*FDPart3_28 - FDPart3_14*FDPart3_46*FDPart3_48*FDPart3_68 - 3.0/2.0*FDPart3_23*FDPart3_7 + FDPart3_31*FDPart3_33*FDPart3_40 - FDPart3_31*FDPart3_33*FDPart3_61*FDPart3_7 + FDPart3_32*FDPart3_33*FDPart3_42 + FDPart3_46*FDPart3_48*FDPart3_54 + FDPart3_47*FDPart3_48*FDPart3_55);
                KDD11[idx] = FDPart3_19*(-7.0/3.0*FDPart3_12*FDPart3_26*puncture_1_P_y + (7.0/6.0)*FDPart3_14*FDPart3_26*puncture_1_P_z - 7.0/3.0*FDPart3_21*FDPart3_4*puncture_0_P_y + (7.0/6.0)*FDPart3_21*FDPart3_7*puncture_0_P_z - FDPart3_22*FDPart3_59 - FDPart3_27*FDPart3_66 + (4.0/3.0)*FDPart3_38 - 2.0/3.0*FDPart3_44 + (4.0/3.0)*FDPart3_52 - 2.0/3.0*FDPart3_57 + FDPart3_60*FDPart3_89 + (4.0/3.0)*FDPart3_62 - 2.0/3.0*FDPart3_64 + FDPart3_67*FDPart3_90 + (4.0/3.0)*FDPart3_69 - 2.0/3.0*FDPart3_71 + FDPart3_97);
                KDD12[idx] = FDPart3_19*(-3.0/2.0*FDPart3_12*FDPart3_26*puncture_1_P_z - FDPart3_12*FDPart3_48*FDPart3_55*FDPart3_68 - FDPart3_14*FDPart3_48*FDPart3_50*FDPart3_68 - 3.0/2.0*FDPart3_14*FDPart3_80 - 3.0/2.0*FDPart3_21*FDPart3_4*puncture_0_P_z + FDPart3_33*FDPart3_34*FDPart3_42 + FDPart3_33*FDPart3_36*FDPart3_40 - FDPart3_33*FDPart3_36*FDPart3_61*FDPart3_7 - FDPart3_33*FDPart3_4*FDPart3_42*FDPart3_61 + FDPart3_48*FDPart3_49*FDPart3_55 + FDPart3_48*FDPart3_50*FDPart3_54 - 3.0/2.0*FDPart3_7*FDPart3_75);
                KDD22[idx] = FDPart3_19*((7.0/6.0)*FDPart3_12*FDPart3_26*puncture_1_P_y - 7.0/3.0*FDPart3_14*FDPart3_26*puncture_1_P_z + (7.0/6.0)*FDPart3_21*FDPart3_4*puncture_0_P_y - 7.0/3.0*FDPart3_21*FDPart3_7*puncture_0_P_z - FDPart3_22*FDPart3_60 - FDPart3_27*FDPart3_67 - 2.0/3.0*FDPart3_38 + (4.0/3.0)*FDPart3_44 - 2.0/3.0*FDPart3_52 + (4.0/3.0)*FDPart3_57 + FDPart3_59*FDPart3_89 - 2.0/3.0*FDPart3_62 + (4.0/3.0)*FDPart3_64 + FDPart3_66*FDPart3_90 - 2.0/3.0*FDPart3_69 + (4.0/3.0)*FDPart3_71 + FDPart3_97);

            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
