//#include "./NRPy_basic_defines.h"
//#include "./NRPy_function_prototypes.h"

// Boundary condition driver routine: Apply BCs to all
// boundary faces of the 3D numerical domain, filling in the
// outer boundary ghost zone layers, starting with the innermost
// layer and working outward.

// Sommerfeld EVOL grid function parameters
const REAL evolgf_at_inf[NUM_EVOL_GFS] = {0.0, 0.0};
const REAL evolgf_radial_falloff_power[NUM_EVOL_GFS] = {3, 3};
REAL evolgf_speed[NUM_EVOL_GFS] = {4.4793122676309114e+06, 4.4793122676309114e+06};


// Function to calculate the radial derivative of a grid function
void conformally_flat_BBH_contraction_term(const paramstruct *restrict params, const int which_gf, const REAL *restrict gfs, REAL *restrict xx[3],
                                           const int8_t FACEXi[3], const int i0, const int i1, const int i2, REAL *restrict _r, REAL *restrict _partial_i_f) {

#include "./conformally_flat_BBH_set_Cparameters.h"

// Initialize derivatives to crazy values, to ensure that
//   we will notice in case they aren't set properly.
REAL fdD0=1e100;
REAL fdD1=1e100;

REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];

int8_t SHIFTSTENCIL0;
int8_t SHIFTSTENCIL1;
int8_t SHIFTSTENCIL2;


// forward/backward finite difference coefficients
const REAL u0 = -49./20.;
const REAL u1 =  6.;
const REAL u2 = -15./2.;
const REAL u3 =  20./3.;
const REAL u4 = -15./4.;
const REAL u5 =  6./5.;
const REAL u6 = -1./6.;

// central finite difference coefficients
const REAL c1 = 3./4.;
const REAL c2 = 3./20.;
const REAL c3 = 1./60;


// On a +x0 or -x0 face, do up/down winding as appropriate:
if(abs(FACEXi[0])==1 || i0+NGHOSTS >= Nxx_plus_2NGHOSTS0 || i0-NGHOSTS <= 0) {
    int8_t SHIFTSTENCIL0 = FACEXi[0];
    if(i0+NGHOSTS >= Nxx_plus_2NGHOSTS0) SHIFTSTENCIL0 = -1;
    if(i0-NGHOSTS <= 0)                  SHIFTSTENCIL0 = +1;
    SHIFTSTENCIL1 = 0;
    SHIFTSTENCIL2 = 0;


    fdD0
        = SHIFTSTENCIL0*(u0*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +u1*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         +u2*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                         +u3*gfs[IDX4S(which_gf,i0+3*SHIFTSTENCIL0,i1+3*SHIFTSTENCIL1,i2+3*SHIFTSTENCIL2)]
                         +u4*gfs[IDX4S(which_gf,i0+4*SHIFTSTENCIL0,i1+4*SHIFTSTENCIL1,i2+4*SHIFTSTENCIL2)]
                         +u5*gfs[IDX4S(which_gf,i0+5*SHIFTSTENCIL0,i1+5*SHIFTSTENCIL1,i2+5*SHIFTSTENCIL2)]
                         +u6*gfs[IDX4S(which_gf,i0+6*SHIFTSTENCIL0,i1+6*SHIFTSTENCIL1,i2+6*SHIFTSTENCIL2)]
                        )*invdx0;

// Not on a +x0 or -x0 face, using centered difference:
} else {
    fdD0 = ( c3*gfs[IDX4S(which_gf,i0+3,i1,i2)]
                         -c2*gfs[IDX4S(which_gf,i0+2,i1,i2)]
                         +c1*gfs[IDX4S(which_gf,i0+1,i1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0-1,i1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0-2,i1,i2)]
                         -c3*gfs[IDX4S(which_gf,i0-3,i1,i2)])*invdx0;
}

// On a +x1 or -x1 face, do up/down winding as appropriate:
if(abs(FACEXi[1])==1 || i1+NGHOSTS >= Nxx_plus_2NGHOSTS1 || i1-NGHOSTS <= 0) {
    int8_t SHIFTSTENCIL1 = FACEXi[1];
    if(i1+NGHOSTS >= Nxx_plus_2NGHOSTS1) SHIFTSTENCIL1 = -1;
    if(i1-NGHOSTS <= 0)                  SHIFTSTENCIL1 = +1;
    SHIFTSTENCIL2 = 0;
    SHIFTSTENCIL0 = 0;


    fdD1
        = SHIFTSTENCIL1*(u0*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +u1*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         +u2*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                         +u3*gfs[IDX4S(which_gf,i0+3*SHIFTSTENCIL0,i1+3*SHIFTSTENCIL1,i2+3*SHIFTSTENCIL2)]
                         +u4*gfs[IDX4S(which_gf,i0+4*SHIFTSTENCIL0,i1+4*SHIFTSTENCIL1,i2+4*SHIFTSTENCIL2)]
                         +u5*gfs[IDX4S(which_gf,i0+5*SHIFTSTENCIL0,i1+5*SHIFTSTENCIL1,i2+5*SHIFTSTENCIL2)]
                         +u6*gfs[IDX4S(which_gf,i0+6*SHIFTSTENCIL0,i1+6*SHIFTSTENCIL1,i2+6*SHIFTSTENCIL2)]
                        )*invdx1;

// Not on a +x1 or -x1 face, using centered difference:
} else {
    fdD1 = ( c3*gfs[IDX4S(which_gf,i0,i1+3,i2)]
                         -c2*gfs[IDX4S(which_gf,i0,i1+2,i2)]
                         +c1*gfs[IDX4S(which_gf,i0,i1+1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0,i1-1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0,i1-2,i2)]
                         -c3*gfs[IDX4S(which_gf,i0,i1-3,i2)])*invdx1;
}

/*
 *  Original SymPy expressions:
 *  "[*_r = sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2),
 *    *_partial_i_f = -fdD0*(-sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)/sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)*cos(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*cos(xx1))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*(-(-sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)/sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)*cos(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*cos(xx1))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))*(AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)) + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)*cos(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 - (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*cos(xx1))*(AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)/(2*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) - AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)))) + fdD1*(AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)/(2*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) - AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*(-(-sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)/sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)*cos(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*cos(xx1))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))*(AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)) + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)*cos(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 - (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*cos(xx1))*(AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)/(2*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(-AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) - AMAX**2*(2*exp(xx0/SINHWAA)/SINHWAA + 2*exp(-xx0/SINHWAA)/SINHWAA)*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*cos(xx1)**2/(2*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2))*cos(xx1)/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)**(3/2))/(sqrt(-(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2/(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2) + 1)*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2))))]"
 */
const double tmp_0 = sin(xx1);
const double tmp_2 = (1.0/(SINHWAA));
const double tmp_4 = exp(tmp_2*xx0);
const double tmp_5 = exp(-tmp_2*xx0);
const double tmp_6 = tmp_4 - tmp_5;
const double tmp_7 = ((AMAX)*(AMAX))/((exp(tmp_2) - exp(-tmp_2))*(exp(tmp_2) - exp(-tmp_2)));
const double tmp_8 = ((tmp_6)*(tmp_6))*tmp_7;
const double tmp_9 = cos(xx1);
const double tmp_11 = ((bScale)*(bScale)) + tmp_8;
const double tmp_12 = tmp_11*((tmp_9)*(tmp_9));
const double tmp_13 = ((tmp_0)*(tmp_0))*tmp_8 + tmp_12;
const double tmp_14 = sqrt(tmp_13);
const double tmp_15 = sqrt(tmp_11);
const double tmp_16 = (1.0/(tmp_14));
const double tmp_18 = -tmp_0*tmp_11*tmp_9 + tmp_0*tmp_8*tmp_9;
const double tmp_19 = pow(tmp_13, -3.0/2.0);
const double tmp_20 = -tmp_0*tmp_15*tmp_16 - tmp_15*tmp_18*tmp_19*tmp_9;
const double tmp_21 = (1.0/sqrt(-tmp_12/tmp_13 + 1));
const double tmp_23 = (1.0/2.0)*tmp_6*tmp_7*(2*tmp_2*tmp_4 + 2*tmp_2*tmp_5);
const double tmp_24 = ((tmp_0)*(tmp_0))*tmp_23 + tmp_23*((tmp_9)*(tmp_9));
const double tmp_25 = -tmp_15*tmp_19*tmp_24*tmp_9 + tmp_16*tmp_23*tmp_9/tmp_15;
const double tmp_26 = (1.0/(tmp_16*tmp_18*tmp_21*tmp_25 - tmp_16*tmp_20*tmp_21*tmp_24));
*_r = tmp_14;
*_partial_i_f = -fdD0*tmp_20*tmp_21*tmp_26 + fdD1*tmp_21*tmp_25*tmp_26;

} // END contraction_term function

void conformally_flat_BBH_apply_bcs_sommerfeld(const paramstruct *restrict params, REAL *restrict xx[3],
                                               const bc_struct *restrict bcstruct, const int NUM_GFS,
                                               const int8_t *restrict gfs_parity, REAL *restrict gfs,
                                               REAL *restrict rhs_gfs) {

    #pragma omp parallel for
        for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
          const REAL char_speed             = evolgf_speed[which_gf];
          const REAL var_at_infinity        = evolgf_at_inf[which_gf];
          const REAL radial_falloff_power = evolgf_radial_falloff_power[which_gf];


#include "./conformally_flat_BBH_set_Cparameters.h"


            for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
                for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
                    const int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
                    const int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
                    const int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;
                    const int8_t FACEX0 = bcstruct->outer[which_gz][pt].FACEi0;
                    const int8_t FACEX1 = bcstruct->outer[which_gz][pt].FACEi1;
                    const int8_t FACEX2 = bcstruct->outer[which_gz][pt].FACEi2;

                    const int8_t FACEXi[3] = {FACEX0, FACEX1, FACEX2};

                    // Initialize derivatives to crazy values, to ensure that
                    //   we will notice in case they aren't set properly.
                    REAL r = 1e100;
                    REAL partial_i_f = 1e100;

                    conformally_flat_BBH_contraction_term(params, which_gf, gfs, xx, FACEXi, i0, i1, i2, &r, &partial_i_f);

                    const REAL invr = 1./r;

                    const REAL source_rhs = -char_speed*(partial_i_f + invr*(gfs[IDX4S(which_gf,i0,i1,i2)] - var_at_infinity));
                    rhs_gfs[IDX4S(which_gf,i0,i1,i2)] = source_rhs;

                    /************* For radial falloff and the extrapolated k term *************/
                    if (radial_falloff_power > 0) {

                      // Move one point away from gz point to compare pure advection to df/dt|interior

                      const int i0_offset = i0+FACEX0;
                      const int i1_offset = i1+FACEX1;
                      const int i2_offset = i2+FACEX2;

                      // Initialize derivatives to crazy values, to ensure that
                      //   we will notice in case they aren't set properly.
                      REAL r_offset = 1e100;
                      REAL partial_i_f_offset = 1e100;

                      conformally_flat_BBH_contraction_term(params, which_gf, gfs, xx, FACEXi, i0_offset, i1_offset, i2_offset, &r_offset, &partial_i_f_offset);

                      const REAL invr_offset = 1./r_offset;

                      // Pure advection: [FIXME: Add equation (appearing in Jupyter notebook documentation)]
                      const REAL extrap_rhs = char_speed*(partial_i_f_offset + invr_offset*(gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] - var_at_infinity));

                      // Take difference between pure advection and df/dt|interior
                      const REAL diff_between_advection_and_f_rhs =
                          rhs_gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] + extrap_rhs;

                      // Solve for k/(r_gz)^n+1 term
                      rhs_gfs[IDX4S(which_gf,i0,i1,i2)] += diff_between_advection_and_f_rhs*pow(r_offset*invr,radial_falloff_power);

                  }
                } // END for(int pt=0;pt<num_ob_gz_pts[which_gz];pt++)

            // Apply INNER (parity) boundary conditions:
                for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
                    const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
                    const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
                    const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
                    const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
                    const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
                    const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;

                    rhs_gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
                            bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * rhs_gfs[IDX4S(which_gf, i0src,i1src,i2src)];
            } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
        } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
} // END function
