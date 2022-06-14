#include "NRPyEllipticET.h"

void NRPyEllipticET_PunctureInitialData_Hyperbolic_Relaxation() {

  DECLARE_CCTK_PARAMETERS;

  // Step 0.a: Check grid parameters
  const int Nxx[3] = { N0,N1,N2 };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }

  // Step 0.b: Print parameter information
  CCTK_INFO("Information about puncture 1:");
  CCTK_VInfo(CCTK_THORNSTRING,"    Bare mass        = %10.3e",puncture_0_bare_mass);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position x       = %10.3e",puncture_0_pos[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position y       = %10.3e",puncture_0_pos[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position z       = %10.3e",puncture_0_pos[2]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum x  = %10.3e",puncture_0_P[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum y  = %10.3e",puncture_0_P[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum z  = %10.3e",puncture_0_P[2]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum x  = %10.3e",puncture_0_S[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum y  = %10.3e",puncture_0_S[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum z  = %10.3e",puncture_0_S[2]);

  CCTK_INFO("Information about puncture 2:");
  CCTK_VInfo(CCTK_THORNSTRING,"    Bare mass        = %10.3e",puncture_1_bare_mass);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position x       = %10.3e",puncture_1_pos[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position y       = %10.3e",puncture_1_pos[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Position z       = %10.3e",puncture_1_pos[2]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum x  = %10.3e",puncture_1_P[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum y  = %10.3e",puncture_1_P[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Lin. momentum z  = %10.3e",puncture_1_P[2]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum x  = %10.3e",puncture_1_S[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum y  = %10.3e",puncture_1_S[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Ang. momentum z  = %10.3e",puncture_1_S[2]);

  // Step 0.c: Set griddata struct
  griddata_struct griddata;
  set_Cparameters_to_default(&griddata.params);

  // Step 0.d: Set parameters from parfile
#include "set_Cparameters_from_parfile.h"

  // Step 0.e: Update parameters if using xy-plane as the orbital plane
  if( CCTK_EQUALS(orbital_plane,"xy") || CCTK_EQUALS(orbital_plane,"yx") ) {
    // The mapping is (x,y,z)_{NRPy} = (y,z,x)_{ETK}.
    // Parameters for puncture #1
    griddata.params.puncture_0_x   = puncture_0_pos[1];
    griddata.params.puncture_0_y   = puncture_0_pos[2];
    griddata.params.puncture_0_z   = puncture_0_pos[0];
    griddata.params.puncture_0_P_x = puncture_0_P[1];
    griddata.params.puncture_0_P_y = puncture_0_P[2];
    griddata.params.puncture_0_P_z = puncture_0_P[0];
    griddata.params.puncture_0_S_x = puncture_0_S[1];
    griddata.params.puncture_0_S_y = puncture_0_S[2];
    griddata.params.puncture_0_S_z = puncture_0_S[0];
    // Parameters for puncture #2
    griddata.params.puncture_1_x   = puncture_1_pos[1];
    griddata.params.puncture_1_y   = puncture_1_pos[2];
    griddata.params.puncture_1_z   = puncture_1_pos[0];
    griddata.params.puncture_1_P_x = puncture_1_P[1];
    griddata.params.puncture_1_P_y = puncture_1_P[2];
    griddata.params.puncture_1_P_z = puncture_1_P[0];
    griddata.params.puncture_1_S_x = puncture_1_S[1];
    griddata.params.puncture_1_S_y = puncture_1_S[2];
    griddata.params.puncture_1_S_z = puncture_1_S[0];
  }

  // Step 0.f: Uniform coordinate grids are stored to *xx[3]
  // Step 0.f.i: Set bcstruct
  {
    int EigenCoord = 1;
    // Step 0.f.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0.g: Find ghostzone mappings; set up bcstruct
    driver_bcstruct(&griddata.params, &griddata.bcstruct, griddata.xx);
    // Step 0.g.i: Free allocated space for xx[][] array
    for(int i=0;i<3;i++) free(griddata.xx[i]);
  }

  // Step 0.h: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
  //           params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
  //           chosen (non-Eigen) CoordSystem.
  int EigenCoord = 0;
  set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);

  // Step 0.i: Time coordinate parameters
  const REAL t_final = domain_size*number_of_LCT;

  // Step 0.j: Set timestep based on smallest proper distance between gridpoints and CFL factor
  REAL dt = find_timestep(&griddata.params, griddata.xx, CFL_FACTOR);

  // Step 0.k: Set maximum number of time steps. Add 0.5 to account
  //          for C rounding down typecasts to integers.
  int N_final = (int)(t_final / dt + 0.5);

  // Step 0.l: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
  //              This is a limitation of the RK method. You are always welcome to declare & allocate
  //              additional gridfunctions by hand.
  if(NUM_AUX_GFS > NUM_EVOL_GFS) {
    printf("Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
    printf("       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
    exit(1);
  }

  // Step 0.m: Declare struct for gridfunctions and allocate memory for y_n_gfs gridfunctions
  MoL_malloc_y_n_gfs(&griddata.params, &griddata.gridfuncs);

  // Step 0.m: Set up precomputed reference metric arrays
  // Step 0.m.i: Allocate space for precomputed reference metric arrays.
  rfm_precompute_rfmstruct_malloc(&griddata.params, &griddata.rfmstruct);

  // Step 0.n.ii: Define precomputed reference metric arrays.
  rfm_precompute_rfmstruct_define(&griddata.params, griddata.xx, &griddata.rfmstruct);

  // Step 1: Set up initial guess
  initial_guess_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.y_n_gfs);

  // Step 2: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  MoL_malloc_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);

  // Step 3: Set up auxevol gridfunctions (psi_background and ADD_times_AUU)
  auxevol_gfs_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.auxevol_gfs);

  // Step 4: Set up wavespeed gridfunction
  wavespeed_gf_all_points(&griddata.params, CFL_FACTOR, dt, griddata.xx, griddata.gridfuncs.auxevol_gfs);

  const int Nxx_plus_2NGHOSTS0 = griddata.params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata.params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata.params.Nxx_plus_2NGHOSTS2;

  // Step 5: Print more useful information
  CCTK_INFO("Grid information:");
  CCTK_VInfo(CCTK_THORNSTRING,"    Points in xx0    =  %d",Nxx[0]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Points in xx1    =  %d",Nxx[1]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Points in xx2    =  %d",Nxx[2]);
  CCTK_VInfo(CCTK_THORNSTRING,"    Domain size      = %10.3e",domain_size);
  CCTK_VInfo(CCTK_THORNSTRING,"    sinh width       = %10.3e",sinh_width);

  CCTK_INFO("Hyperbolic relaxation information:");
  CCTK_VInfo(CCTK_THORNSTRING,"    Damping strength = %10.3e",eta_damping);
  CCTK_VInfo(CCTK_THORNSTRING,"    Target residual  = %10.3e",pow(10.0,target_log_residual));
  CCTK_VInfo(CCTK_THORNSTRING,"    Time step        = %10.3e",dt);
  CCTK_VInfo(CCTK_THORNSTRING,"    Final time       = %10.3e",number_of_LCT*domain_size);
  CCTK_VInfo(CCTK_THORNSTRING,"    Max iterations   =  %d"   ,N_final);

  // Step 6: Relaxation time loop
  // Step 6.a: Start timing for NRPyEllipticET
  struct timespec start, end;
  clock_gettime(CLOCK_REALTIME, &start);

  // Step 6.b: Set frequency with which the L2-norm of residual is computed
  int output_l2_norm_residual_every_N = (int)((REAL)N_final/1000.0);
  if( output_l2_norm_residual_every_N == 0 ) output_l2_norm_residual_every_N = 1;

  for(int n=0;n<=N_final;n++) { // Main loop to progress forward in time.

    // Step 6.c: Set current time to correct value
    griddata.params.time = ((REAL)n)*dt;

    // Step 6.d: Step forward one timestep (t -> t+dt) in time using
    //           chosen RK-like MoL timestepping algorithm
    MoL_step_forward_in_time(&griddata, dt);

    // Step 6.e: Check convergence & print information
    if( (n%output_l2_norm_residual_every_N == 0) || (n==N_final) ) {

      /*
      // Step 6.e.i: Compute residual and store it at diagnostic_output_gfs(UUFG)
      residual_all_points(&griddata.params, &griddata.rfmstruct, griddata.gridfuncs.auxevol_gfs,
                          griddata.gridfuncs.y_n_gfs, griddata.gridfuncs.diagnostic_output_gfs);

      const int residual_gf_index   = UUGF;
      const REAL integration_radius = 100.0;
      const REAL l2_norm_residual   = L2_norm_of_gf(&griddata.params, residual_gf_index, integration_radius, griddata.xx,
                                                    griddata.gridfuncs.diagnostic_output_gfs);
      */

      // Step 6.e.ii: Compute time elapsed during the executation of this thorn
      clock_gettime(CLOCK_REALTIME, &end);
      long long unsigned int time_elapsed_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
      REAL time_elapsed_s = (REAL)time_elapsed_ns*1e-9;

      // Step 6.e.iii: Compute percentage of the run that has completed
      REAL completion = griddata.params.time/t_final;

      // Step 6.e.iv: Estimate how much time is left
      REAL time_left_s = time_elapsed_s*(1.0/completion - 1.0);

      // Step 6.e.v: Print information to the user
      printf("INFO (NRPyEllipticET): It: %05d, t: %4.2lf, completion: %5.1lf%%, runtime: %3.0lf s, ETA: %3.0lf s\r",
             n,griddata.params.time,completion*100,time_elapsed_s,time_left_s);
    }
  } // End main loop to progress forward in time.

  // Step 7: Initialize NRPyEllipticET_xx
  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) NRPyEllipticET_xx[0][i0] = griddata.xx[0][i0];
  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) NRPyEllipticET_xx[1][i1] = griddata.xx[1][i1];
  for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++) NRPyEllipticET_xx[2][i2] = griddata.xx[2][i2];

  // Step 8: Initialize NRPyEllipticET_uu
#pragma omp parallel for
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0,
              0,Nxx_plus_2NGHOSTS1,
              0,Nxx_plus_2NGHOSTS2) {
    NRPyEllipticET_uu[IDX4S(UUGF,i0,i1,i2)] = griddata.gridfuncs.y_n_gfs[IDX4S(UUGF,i0,i1,i2)];
  }

  // Step 9: Print final information about the run
  clock_gettime(CLOCK_REALTIME, &end);
  const long long unsigned int time_in_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
  CCTK_VInfo(CCTK_THORNSTRING,"Total execution time: %.0lf s", ((REAL) time_in_ns )*1.0e-9);

  // Step 10: Free all allocated memory
  rfm_precompute_rfmstruct_freemem(&griddata.params, &griddata.rfmstruct);
  freemem_bcstruct(&griddata.params, &griddata.bcstruct);
  MoL_free_memory_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  MoL_free_memory_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
}
