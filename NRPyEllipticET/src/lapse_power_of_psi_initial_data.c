#include "NRPyEllipticET.h"

/*
 * Lapse initial data: power of psi
 */
void lapse_power_of_psi_initial_data(const cGH *restrict cctkGH,
                                     const paramstruct *restrict params,
                                     const REAL *restrict x,
                                     const REAL *restrict y,
                                     const REAL *restrict z,
                                     const REAL *restrict in_gfs,
                                     REAL *restrict alpha) {
#include "./set_Cparameters.h"

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
                const double FDPart3_1 = -Carty;
                alpha[idx] = pow((1.0/2.0)*puncture_0_bare_mass/sqrt(((-FDPart3_1 - puncture_0_y)*(-FDPart3_1 - puncture_0_y)) + ((-puncture_0_x + Cartx)*(-puncture_0_x + Cartx)) + ((-puncture_0_z + Cartz)*(-puncture_0_z + Cartz))) + (1.0/2.0)*puncture_1_bare_mass/sqrt(((-FDPart3_1 - puncture_1_y)*(-FDPart3_1 - puncture_1_y)) + ((-puncture_1_x + Cartx)*(-puncture_1_x + Cartx)) + ((-puncture_1_z + Cartz)*(-puncture_1_z + Cartz))) + uu + 1, lapse_exponent);

            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
