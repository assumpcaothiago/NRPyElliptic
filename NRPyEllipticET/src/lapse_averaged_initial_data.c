#include "NRPyEllipticET.h"

/*
 * Lapse initial data: averaged
 */
void lapse_averaged_initial_data(const cGH *restrict cctkGH,
                                 const char *restrict orbital_plane,
                                 const paramstruct *restrict params,
                                 const REAL *restrict x,
                                 const REAL *restrict y,
                                 const REAL *restrict z,
                                 const REAL *restrict position_shift,
                                 const REAL *restrict in_gfs,
                                 REAL *restrict alpha) {
#include "./set_Cparameters.h"

    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {

                const int idx      = CCTK_GFINDEX3D(cctkGH,i0,i1,i2);

                // Set the (x,y,z) Cartesian coordinates
                double Cartx = x[idx];
                double Carty = y[idx];
                double Cartz = z[idx];

                /*
                 * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
                 */
                const double FDPart3_1 = -Carty;
                const double FDPart3_3 = (1.0/2.0)*puncture_0_bare_mass/sqrt(((-FDPart3_1 - puncture_0_y)*(-FDPart3_1 - puncture_0_y)) + ((-puncture_0_x + Cartx)*(-puncture_0_x + Cartx)) + ((-puncture_0_z + Cartz)*(-puncture_0_z + Cartz))) + (1.0/2.0)*puncture_1_bare_mass/sqrt(((-FDPart3_1 - puncture_1_y)*(-FDPart3_1 - puncture_1_y)) + ((-puncture_1_x + Cartx)*(-puncture_1_x + Cartx)) + ((-puncture_1_z + Cartz)*(-puncture_1_z + Cartz)));
                alpha[idx] = (1.0/2.0)*(1 - FDPart3_3)/(FDPart3_3 + 1) + 1.0/2.0;

            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
