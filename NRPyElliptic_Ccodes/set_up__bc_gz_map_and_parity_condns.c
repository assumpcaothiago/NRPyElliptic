#include "././NRPy_basic_defines.h"
#include "././NRPy_function_prototypes.h"
void set_up__bc_gz_map_and_parity_condns(const paramstruct *restrict params,
                                             REAL *restrict xx[3], gz_map *restrict bc_gz_map,parity_condition *restrict bc_parity_conditions) {
#include "./set_Cparameters.h"


  // xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx0;
  // -> xxmin[0] = xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0
  const REAL xxmin[3] = { xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
                          xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
                          xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };
  //fprintf(stderr,"hey inside setbc: %e %e %e | %e %e\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],dxx0);
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0,0,Nxx_plus_2NGHOSTS1,0,Nxx_plus_2NGHOSTS2) {
    // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
    REAL xCart[3];
    {
      // xx_to_Cart for EigenCoordinate SymTP (orig coord = SinhSymTP):
      REAL xx0 = xx[0][i0];
      REAL xx1 = xx[1][i1];
      REAL xx2 = xx[2][i2];
      /*
       *  Original SymPy expressions:
       *  "[xCart[0] = xx0*sin(xx1)*cos(xx2),
       *    xCart[1] = xx0*sin(xx1)*sin(xx2),
       *    xCart[2] = sqrt(bScale**2 + xx0**2)*cos(xx1)]"
       */
      {
        const double tmp_0 = xx0*sin(xx1);
        xCart[0] = tmp_0*cos(xx2);
        xCart[1] = tmp_0*sin(xx2);
        xCart[2] = sqrt(((bScale)*(bScale)) + ((xx0)*(xx0)))*cos(xx1);
      }
    }

    //EigenCoord_xx_to_Cart(params, xx, i0,i1,i2, xCart);
    REAL Cartx = xCart[0];
    REAL Carty = xCart[1];
    REAL Cartz = xCart[2];

    // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
    //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
    //      the point is an outer boundary point.
    //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
    //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
    //      appropriate parity condition (+/- 1).
    REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
    // Cart_to_xx for EigenCoordinate SymTP (orig coord = SinhSymTP);
    /*
     *  Original SymPy expressions:
     *  "[Cart_to_xx0_inbounds = M_SQRT1_2*sqrt(Cartx**2 + Carty**2 + Cartz**2 - bScale**2 + sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)),
     *    Cart_to_xx1_inbounds = acos(M_SQRT1_2*sqrt(1 + (Cartx**2 + Carty**2 + Cartz**2)/bScale**2 - sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)/bScale**2)*sign(Cartz)),
     *    Cart_to_xx2_inbounds = atan2(Carty, Cartx)]"
     */
    {
      const double tmp_1 = ((bScale)*(bScale));
      const double tmp_2 = ((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz));
      const double tmp_3 = sqrt(-4*((Cartz)*(Cartz))*tmp_1 + ((bScale)*(bScale)*(bScale)*(bScale)) + 2*tmp_1*tmp_2 + ((tmp_2)*(tmp_2)));
      const double tmp_4 = (1.0/(tmp_1));
      Cart_to_xx0_inbounds = M_SQRT1_2*sqrt(-tmp_1 + tmp_2 + tmp_3);
      Cart_to_xx1_inbounds = acos(M_SQRT1_2*sqrt(tmp_2*tmp_4 - tmp_3*tmp_4 + 1)*(((Cartz) > 0) - ((Cartz) < 0)));
      Cart_to_xx2_inbounds = atan2(Carty, Cartx);
    }

    const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
    const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
    const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );

    // Step 2.a: (Sanity/validation check) Convert the interior point
    //           x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) to Cartesian coordinates,
    //           make sure that the Cartesian coordinate matches the Cartesian coordinate of
    //           x0(i0),x1(i1),x2(i2). If not, error out!
    REAL xCart_orig[3]; for(int ii=0;ii<3;ii++) xCart_orig[ii] = xCart[ii];
    //EigenCoord_xx_to_Cart(params, xx, i0_inbounds,i1_inbounds,i2_inbounds, xCart);
    {
      // xx_to_Cart for EigenCoordinate SymTP (orig coord = SinhSymTP):
      REAL xx0 = xx[0][i0_inbounds];
      REAL xx1 = xx[1][i1_inbounds];
      REAL xx2 = xx[2][i2_inbounds];
      /*
       *  Original SymPy expressions:
       *  "[xCart[0] = xx0*sin(xx1)*cos(xx2),
       *    xCart[1] = xx0*sin(xx1)*sin(xx2),
       *    xCart[2] = sqrt(bScale**2 + xx0**2)*cos(xx1)]"
       */
      {
        const double tmp_0 = xx0*sin(xx1);
        xCart[0] = tmp_0*cos(xx2);
        xCart[1] = tmp_0*sin(xx2);
        xCart[2] = sqrt(((bScale)*(bScale)) + ((xx0)*(xx0)))*cos(xx1);
      }
    }


//fprintf(stderr,"Cartesian agreement: ( %.15e %.15e %.15e ) ?= ( %.15e %.15e %.15e )\n",
// (double)xCart_orig[0],(double)xCart_orig[1],(double)xCart_orig[2],
// (double)xCart[0],(double)xCart[1],(double)xCart[2]);

#define EPS_ABS 1e-8
    if(fabs( (double)(xCart_orig[0] - xCart[0]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[1] - xCart[1]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[2] - xCart[2]) ) > EPS_ABS) {
       fprintf(stderr,"Error. SinhSymTP: Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
               (double)xCart_orig[0],(double)xCart_orig[1],(double)xCart_orig[2],
               (double)xCart[0],(double)xCart[1],(double)xCart[2],
               xx[0][i0],xx[1][i1],xx[2][i2],
               xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
               Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
      exit(1);
    }

    // Step 3: Set bc_gz_map and bc_parity_conditions.
    if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0) {
      // Step 3.a: Iff we are on an outer boundary point or in the grid
      //           interior, i0_inbounds==i0, i1_inbounds==i1, and
      //           i2_inbounds==i2, and inner boundary conditions do not
      //           apply: set bc_gz_map to -1, and parity=1.
      bc_gz_map[IDX3S(i0,i1,i2)].i0=-1;
      bc_gz_map[IDX3S(i0,i1,i2)].i1=-1;
      bc_gz_map[IDX3S(i0,i1,i2)].i2=-1;
      for(int which_parity=0; which_parity<10; which_parity++) {
        bc_parity_conditions[IDX3S(i0,i1,i2)].parity[which_parity] = 1;
      }
    } else {
      // Step 3.b: If we are on an *inner* boundary point:
      // 1. Set bc_gz_map at (i0,i1,i2) to the point
      //    in the interior to which this boundary
      //    point maps, and
      // 2. Perform the unit vector dot products
      //    necessary to set all 10 possible parity
      //    conditions, calling function
      //    set_parity_from_unit_vector_dot_product()
      bc_gz_map[IDX3S(i0,i1,i2)].i0=i0_inbounds;
      bc_gz_map[IDX3S(i0,i1,i2)].i1=i1_inbounds;
      bc_gz_map[IDX3S(i0,i1,i2)].i2=i2_inbounds;
      const REAL xx0 = xx[0][i0];
      const REAL xx1 = xx[1][i1];
      const REAL xx2 = xx[2][i2];
      const REAL xx0_inbounds = xx[0][i0_inbounds];
      const REAL xx1_inbounds = xx[1][i1_inbounds];
      const REAL xx2_inbounds = xx[2][i2_inbounds];
      REAL REAL_parity_array[10];
      {
      // Evaluate dot products needed for setting parity
      //     conditions at a given point (xx0,xx1,xx2),
      //     using C code generated by NRPy+

      // NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
      //      ten parity conditions, in given coordinate system.
      //      Needed for automatically determining sign of tensor across coordinate boundary.
      // Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
        /*
         *  Original SymPy expressions:
         *  "[REAL_parity_array[0] = 1,
         *    REAL_parity_array[1] = AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)),
         *    REAL_parity_array[2] = AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)),
         *    REAL_parity_array[3] = sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds),
         *    REAL_parity_array[4] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))**2,
         *    REAL_parity_array[5] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
         *    REAL_parity_array[6] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
         *    REAL_parity_array[7] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))**2,
         *    REAL_parity_array[8] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
         *    REAL_parity_array[9] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))**2]"
         */
        {
          const double tmp_0 = (1.0/(SINHWAA));
          const double tmp_2 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
          const double tmp_4 = exp(tmp_0*xx0_inbounds) - exp(-tmp_0*xx0_inbounds);
          const double tmp_5 = ((AMAX)*(AMAX))/((exp(tmp_0) - exp(-tmp_0))*(exp(tmp_0) - exp(-tmp_0)));
          const double tmp_6 = ((bScale)*(bScale));
          const double tmp_7 = sin(xx1);
          const double tmp_8 = ((tmp_2)*(tmp_2))*tmp_5;
          const double tmp_9 = sin(xx1_inbounds);
          const double tmp_10 = ((tmp_4)*(tmp_4))*tmp_5;
          const double tmp_11 = 1/(sqrt(tmp_10 + tmp_6*((tmp_9)*(tmp_9)))*sqrt(tmp_6*((tmp_7)*(tmp_7)) + tmp_8));
          const double tmp_12 = tmp_11*tmp_2*tmp_4*tmp_5*cos(xx1)*cos(xx1_inbounds);
          const double tmp_13 = sin(xx2)*sin(xx2_inbounds);
          const double tmp_14 = tmp_11*tmp_7*tmp_9*sqrt(tmp_10 + tmp_6)*sqrt(tmp_6 + tmp_8);
          const double tmp_15 = cos(xx2)*cos(xx2_inbounds);
          const double tmp_16 = tmp_12 + tmp_13*tmp_14 + tmp_14*tmp_15;
          const double tmp_17 = tmp_12*tmp_13 + tmp_12*tmp_15 + tmp_14;
          const double tmp_18 = tmp_13 + tmp_15;
          REAL_parity_array[0] = 1;
          REAL_parity_array[1] = tmp_16;
          REAL_parity_array[2] = tmp_17;
          REAL_parity_array[3] = tmp_18;
          REAL_parity_array[4] = ((tmp_16)*(tmp_16));
          REAL_parity_array[5] = tmp_16*tmp_17;
          REAL_parity_array[6] = tmp_16*tmp_18;
          REAL_parity_array[7] = ((tmp_17)*(tmp_17));
          REAL_parity_array[8] = tmp_17*tmp_18;
          REAL_parity_array[9] = ((tmp_18)*(tmp_18));
        }

      }
      //eval_symbolic_dot_products_to_set_parity_conditions(params,  REAL_parity_array,  xx0,xx1,xx2,
      //                                                          xx0_inbounds,xx1_inbounds,xx2_inbounds);
      for(int whichparity=0;whichparity<10;whichparity++) {
          //printf("Good? Parity %d evaluated to %e\n",whichparity,(double)REAL_parity_array[whichparity]);
          // Perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
          if( (REAL_parity_array[whichparity]  > 0 && fabs(REAL_parity_array[whichparity] - (+1)) > 1e-8) ||
              (REAL_parity_array[whichparity] <= 0 && fabs(REAL_parity_array[whichparity] - (-1)) > 1e-8) ) {
              fprintf(stderr,"Error at point (%d %d %d); (%e %e %e); maps to (%e %e %e).\n",
                      i0,i1,i2, xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
              fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
                      REAL_parity_array[whichparity]);
              exit(1);
          }
          if(REAL_parity_array[whichparity] < 0.0) bc_parity_conditions[IDX3S(i0,i1,i2)].parity[whichparity] = -1;
          if(REAL_parity_array[whichparity] > 0.0) bc_parity_conditions[IDX3S(i0,i1,i2)].parity[whichparity] = +1;
      } // END for(int whichparity=0;whichparity<10;whichparity++)
    } // END if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0)
  } // END LOOP_REGION(0,Nxx_plus_2NGHOSTS0,0,Nxx_plus_2NGHOSTS1,0,Nxx_plus_2NGHOSTS2)
}
