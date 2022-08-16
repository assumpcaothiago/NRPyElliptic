#include "../NRPyEllipticET.h"
#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"

void NRPyEllipticET_conformally_flat_BBH_Initialize_ADMBase( const cGH *restrict cctkGH,
                                                             const int *restrict cctk_lsh,
                                                             const CCTK_REAL *restrict x,
                                                             const CCTK_REAL *restrict y,
                                                             const CCTK_REAL *restrict z,
                                                             CCTK_REAL *restrict uuGF,
                                                             CCTK_REAL *restrict alp,
                                                             CCTK_REAL *restrict betax,
                                                             CCTK_REAL *restrict betay,
                                                             CCTK_REAL *restrict betaz,
                                                             CCTK_REAL *restrict gxx,
                                                             CCTK_REAL *restrict gxy,
                                                             CCTK_REAL *restrict gxz,
                                                             CCTK_REAL *restrict gyy,
                                                             CCTK_REAL *restrict gyz,
                                                             CCTK_REAL *restrict gzz,
                                                             CCTK_REAL *restrict kxx,
                                                             CCTK_REAL *restrict kxy,
                                                             CCTK_REAL *restrict kxz,
                                                             CCTK_REAL *restrict kyy,
                                                             CCTK_REAL *restrict kyz,
                                                             CCTK_REAL *restrict kzz ) {

  // Step 0: Declare CCTK arguments (which includes the gridfunctions) and parameters
  DECLARE_CCTK_PARAMETERS;

  paramstruct params;
#include "conformally_flat_BBH_set_Cparameters_from_parfile_paramstruct.h"

  // Shift position of the punctures
  params.puncture_0_x += position_shift[0];
  params.puncture_0_y += position_shift[1];
  params.puncture_0_z += position_shift[2];
  params.puncture_1_x += position_shift[0];
  params.puncture_1_y += position_shift[1];
  params.puncture_1_z += position_shift[2];

  const int Nxx_plus_2NGHOSTS0    = N0 + 2*NGHOSTS;
  const int Nxx_plus_2NGHOSTS1    = N1 + 2*NGHOSTS;
  const int Nxx_plus_2NGHOSTS2    = N2 + 2*NGHOSTS;

  // Step 2: Interpolate data to the ETK grid
  // Step 2b: Set the origin and step sizes of the grids
  const CCTK_REAL origin[3] = { NRPyEllipticET_xx[0][0],NRPyEllipticET_xx[1][0],NRPyEllipticET_xx[2][0] };
  const CCTK_REAL deltas[3] = {
    NRPyEllipticET_xx[0][1]-NRPyEllipticET_xx[0][0],
    NRPyEllipticET_xx[1][1]-NRPyEllipticET_xx[1][0],
    NRPyEllipticET_xx[2][1]-NRPyEllipticET_xx[2][0]
  };

  const CCTK_INT input_array_dims[3] = { Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2 };

  // Step 3: Loop over the CCTK grid and interpolate the gridfunctions
  int imin = 0; int imax = cctk_lsh[0];
  int jmin = 0; int jmax = cctk_lsh[1];
  int kmin = 0; int kmax = cctk_lsh[2];

  // Set (xx0,xx1,xx2)
  const int num_pts = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  CCTK_REAL *points_dest_grid_xx0 = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_pts);
  CCTK_REAL *points_dest_grid_xx1 = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_pts);
  CCTK_REAL *points_dest_grid_xx2 = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_pts);

  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        // Local index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Read in the points (x,y,z) of the Cartesian grid
        CCTK_REAL xCart[3];
        if( CCTK_EQUALS(orbital_plane,"xy") || CCTK_EQUALS(orbital_plane,"yx") ) {
          xCart[0] = y[index] - position_shift[1];
          xCart[1] = z[index] - position_shift[2];
          xCart[2] = x[index] - position_shift[0];
        }
        else if( CCTK_EQUALS(orbital_plane,"xz") || CCTK_EQUALS(orbital_plane,"zx") ) {
          xCart[0] = x[index] - position_shift[0];
          xCart[1] = y[index] - position_shift[1];
          xCart[2] = z[index] - position_shift[2];
        }
        else {
          CCTK_ERROR("Unsupported orbital plane. Available options are \"xy\" (\"yx\") or \"xz\" (\"zx\").");
        }
        // Now get (xx0,xx1,xx2) from (x,y,z)
        CCTK_REAL xx_star[3];
        conformally_flat_BBH_get_xx_from_xyz(domain_size,sinh_width,foci_position,xCart,xx_star);

        if( CCTK_isnan(xx_star[0]*xx_star[1]*xx_star[2]) ) {
          // Found a NAN: probably too close to the puncture. Let's try again.
          xCart[0] += NRPy_epsilon;
          xCart[1] += NRPy_epsilon;
          xCart[2] += NRPy_epsilon;

          conformally_flat_BBH_get_xx_from_xyz(domain_size,sinh_width,foci_position,xCart,xx_star);

          if( CCTK_isnan(xx_star[0]*xx_star[1]*xx_star[2]) ) {
            // Still found a NAN. Error out.
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "NAN found: %e %e %e from %e %e %e. Please considering increasing the value of NRPy_epsilon.\n",
                        xx_star[0],xx_star[1],xx_star[2],xCart[0],xCart[1],xCart[2]);
          }
        }

        points_dest_grid_xx0[index] = xx_star[0];
        points_dest_grid_xx1[index] = xx_star[1];
        points_dest_grid_xx2[index] = xx_star[2];

        if( CCTK_isnan(xx_star[0]*xx_star[1]*xx_star[2]) ) {
          printf("NAN found: %e %e %e from %e %e %e\n",xx_star[0],xx_star[1],xx_star[2],xCart[0],xCart[1],xCart[2]);
        }

      } // for(int i=imin;i<imax;i++)
    } // for(int j=jmin;j<jmax;j++)
  } // for(int k=kmin;k<kmax;k++)

  // Set pointer to gridpoints
  const void* interp_coords[3] = {
    (const void *) points_dest_grid_xx0,
    (const void *) points_dest_grid_xx1,
    (const void *) points_dest_grid_xx2
  };

  // Input gridfunctions
  CCTK_REAL input_array[Nxx_plus_2NGHOSTS2][Nxx_plus_2NGHOSTS1][Nxx_plus_2NGHOSTS0];
#pragma omp parallel for
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0,
              0,Nxx_plus_2NGHOSTS1,
              0,Nxx_plus_2NGHOSTS2) {
    input_array[i2][i1][i0] = NRPyEllipticET_uu[IDX4S(UUGF,i0,i1,i2)];
  }

  // Output gridfunctions
  CCTK_REAL* output_array = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_pts);

  CCTK_INFO("Interpolating puncture initial data");

  NRPyEllipticET_conformally_flat_BBH_interpolate_solution_to_ADMBase(num_pts,
                                                                      input_array_dims,
                                                                      origin,deltas,
                                                                      interp_coords,
                                                                      input_array,
                                                                      output_array);

  // Change the params to conform to the ETK grid
  params.Nxx_plus_2NGHOSTS0 = cctk_lsh[0];
  params.Nxx_plus_2NGHOSTS1 = cctk_lsh[1];
  params.Nxx_plus_2NGHOSTS2 = cctk_lsh[2];

  CCTK_INFO("Computing ADMBase quantities from NRPyEllipticET_uu_gf");

#pragma omp parallel for
  LOOP_REGION(0,cctk_lsh[0],
              0,cctk_lsh[1],
              0,cctk_lsh[2]) {
    uuGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = output_array[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)];
  }

  // Compute ADM quantities from uu
  conformally_flat_BBH_compute_ADM_Cartesian_quantities_from_uu(cctkGH,orbital_plane,&params,
                                                                x,y,z,position_shift,
                                                                output_array,
                                                                alp,betax,betay,betaz,
                                                                gxx,gxy,gxz,gyy,gyz,gzz,
                                                                kxx,kxy,kxz,kyy,kyz,kzz);

  // Free memory
  free(points_dest_grid_xx0);
  free(points_dest_grid_xx1);
  free(points_dest_grid_xx2);
  free((void *)output_array);

}
