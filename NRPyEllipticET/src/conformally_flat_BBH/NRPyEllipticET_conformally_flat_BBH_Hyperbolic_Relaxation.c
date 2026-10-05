#include "../NRPyEllipticET.h"
#include "./conformally_flat_BBH_NRPy_basic_defines.h"
#include "./conformally_flat_BBH_NRPy_function_prototypes.h"

void NRPyEllipticET_conformally_flat_BBH_Hyperbolic_Relaxation() {

  DECLARE_CCTK_PARAMETERS;

  // Set the size of the local grid
  const CCTK_INT Nxx_plus_2NGHOSTS0    = N0 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS1    = N1 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS2    = N2 + 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  // Allocate memory for NRPyEllipticET_xx and NRPyEllipticET_uu
  NRPyEllipticET_xx[0] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS0);
  NRPyEllipticET_xx[1] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS1);
  NRPyEllipticET_xx[2] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS2);
  NRPyEllipticET_uu    = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*Nxx_plus_2NGHOSTS_tot);

  // Check everything is okay
  if( !NRPyEllipticET_uu || !NRPyEllipticET_xx[0] || !NRPyEllipticET_xx[1] || !NRPyEllipticET_xx[2] ) {
    CCTK_ERROR("Could not allocate memory for NRPyEllipticET_uu or NRPyEllipticET_xx!");
  }

  // Step 0.a: Check grid parameters
  const int Nxx[3] = { N0,N1,N2 };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    CCTK_ERROR("Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number."
               " For example, in case of angular directions, proper symmetry zones will not exist.");
  }

  // Step 0.b: Print parameter information
  CCTK_INFO("Information about puncture 1:");
  CCTK_VINFO("    Bare mass        = %10.3e",conformally_flat_BBH_puncture_0_bare_mass);
  CCTK_VINFO("    Position x       = %10.3e",conformally_flat_BBH_puncture_0_pos[0]);
  CCTK_VINFO("    Position y       = %10.3e",conformally_flat_BBH_puncture_0_pos[1]);
  CCTK_VINFO("    Position z       = %10.3e",conformally_flat_BBH_puncture_0_pos[2]);
  CCTK_VINFO("    Lin. momentum x  = %10.3e",conformally_flat_BBH_puncture_0_P[0]);
  CCTK_VINFO("    Lin. momentum y  = %10.3e",conformally_flat_BBH_puncture_0_P[1]);
  CCTK_VINFO("    Lin. momentum z  = %10.3e",conformally_flat_BBH_puncture_0_P[2]);
  CCTK_VINFO("    Ang. momentum x  = %10.3e",conformally_flat_BBH_puncture_0_S[0]);
  CCTK_VINFO("    Ang. momentum y  = %10.3e",conformally_flat_BBH_puncture_0_S[1]);
  CCTK_VINFO("    Ang. momentum z  = %10.3e",conformally_flat_BBH_puncture_0_S[2]);

  CCTK_INFO("Information about puncture 2:");
  CCTK_VINFO("    Bare mass        = %10.3e",conformally_flat_BBH_puncture_1_bare_mass);
  CCTK_VINFO("    Position x       = %10.3e",conformally_flat_BBH_puncture_1_pos[0]);
  CCTK_VINFO("    Position y       = %10.3e",conformally_flat_BBH_puncture_1_pos[1]);
  CCTK_VINFO("    Position z       = %10.3e",conformally_flat_BBH_puncture_1_pos[2]);
  CCTK_VINFO("    Lin. momentum x  = %10.3e",conformally_flat_BBH_puncture_1_P[0]);
  CCTK_VINFO("    Lin. momentum y  = %10.3e",conformally_flat_BBH_puncture_1_P[1]);
  CCTK_VINFO("    Lin. momentum z  = %10.3e",conformally_flat_BBH_puncture_1_P[2]);
  CCTK_VINFO("    Ang. momentum x  = %10.3e",conformally_flat_BBH_puncture_1_S[0]);
  CCTK_VINFO("    Ang. momentum y  = %10.3e",conformally_flat_BBH_puncture_1_S[1]);
  CCTK_VINFO("    Ang. momentum z  = %10.3e",conformally_flat_BBH_puncture_1_S[2]);

  // Step 0.c: Set griddata struct
  griddata_struct griddata;
  conformally_flat_BBH_set_Cparameters_to_default(&griddata.params);

  // Step 0.d: Set parameters from parfile
#include "conformally_flat_BBH_set_Cparameters_from_parfile.h"

  // Step 0.e: Update parameters if using xy-plane as the orbital plane
  if( CCTK_EQUALS(orbital_plane,"xy") || CCTK_EQUALS(orbital_plane,"yx") ) {
    // The mapping is (x,y,z)_{NRPy} = (y,z,x)_{ETK}.
    // Parameters for puncture #1
    griddata.params.puncture_0_x   = conformally_flat_BBH_puncture_0_pos[1];
    griddata.params.puncture_0_y   = conformally_flat_BBH_puncture_0_pos[2];
    griddata.params.puncture_0_z   = conformally_flat_BBH_puncture_0_pos[0];
    griddata.params.puncture_0_P_x = conformally_flat_BBH_puncture_0_P[1];
    griddata.params.puncture_0_P_y = conformally_flat_BBH_puncture_0_P[2];
    griddata.params.puncture_0_P_z = conformally_flat_BBH_puncture_0_P[0];
    griddata.params.puncture_0_S_x = conformally_flat_BBH_puncture_0_S[1];
    griddata.params.puncture_0_S_y = conformally_flat_BBH_puncture_0_S[2];
    griddata.params.puncture_0_S_z = conformally_flat_BBH_puncture_0_S[0];
    // Parameters for puncture #2
    griddata.params.puncture_1_x   = conformally_flat_BBH_puncture_1_pos[1];
    griddata.params.puncture_1_y   = conformally_flat_BBH_puncture_1_pos[2];
    griddata.params.puncture_1_z   = conformally_flat_BBH_puncture_1_pos[0];
    griddata.params.puncture_1_P_x = conformally_flat_BBH_puncture_1_P[1];
    griddata.params.puncture_1_P_y = conformally_flat_BBH_puncture_1_P[2];
    griddata.params.puncture_1_P_z = conformally_flat_BBH_puncture_1_P[0];
    griddata.params.puncture_1_S_x = conformally_flat_BBH_puncture_1_S[1];
    griddata.params.puncture_1_S_y = conformally_flat_BBH_puncture_1_S[2];
    griddata.params.puncture_1_S_z = conformally_flat_BBH_puncture_1_S[0];
  }

  // Step 0.f: Uniform coordinate grids are stored to *xx[3]
  // Step 0.f.i: Set bcstruct
  {
    int EigenCoord = 1;
    // Step 0.f.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    conformally_flat_BBH_set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0.g: Find ghostzone mappings; set up bcstruct
    conformally_flat_BBH_driver_bcstruct(&griddata.params, &griddata.bcstruct, griddata.xx);
    // Step 0.g.i: Free allocated space for xx[][] array
    for(int i=0;i<3;i++) free(griddata.xx[i]);
  }

  // Step 0.h: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
  //           params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
  //           chosen (non-Eigen) CoordSystem.
  int EigenCoord = 0;
  conformally_flat_BBH_set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);

  // Step 0.i: Set timestep based on smallest proper distance between gridpoints and CFL factor
  REAL dt = conformally_flat_BBH_find_timestep(&griddata.params, griddata.xx, CFL_FACTOR);

  // Step 0.j: Set maximum number of time steps
  int N_final = max_iterations;

  // Step 0.l: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
  //              This is a limitation of the RK method. You are always welcome to declare & allocate
  //              additional gridfunctions by hand.
  if(NUM_AUX_GFS > NUM_EVOL_GFS) {
    CCTK_ERROR("NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,"
               " or allocate (malloc) by hand storage for *diagnostic_output_gfs.");
  }

  // Step 0.m: Declare struct for gridfunctions and allocate memory for y_n_gfs gridfunctions
  conformally_flat_BBH_MoL_malloc_y_n_gfs(&griddata.params, &griddata.gridfuncs);

  // Step 0.m: Set up precomputed reference metric arrays
  // Step 0.m.i: Allocate space for precomputed reference metric arrays.
  conformally_flat_BBH_rfm_precompute_rfmstruct_malloc(&griddata.params, &griddata.rfmstruct);

  // Step 0.n.ii: Define precomputed reference metric arrays.
  conformally_flat_BBH_rfm_precompute_rfmstruct_define(&griddata.params, griddata.xx, &griddata.rfmstruct);

  // Step 1: Set up initial guess
  conformally_flat_BBH_initial_guess_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.y_n_gfs);

  // Step 2: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  conformally_flat_BBH_MoL_malloc_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);

  // Step 3: Set up auxevol gridfunctions (psi_background and ADD_times_AUU)
  conformally_flat_BBH_auxevol_gfs_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.auxevol_gfs);

  // Step 4: Set up wavespeed gridfunction
  conformally_flat_BBH_wavespeed_gf_all_points(&griddata.params, CFL_FACTOR, dt, griddata.xx, griddata.gridfuncs.auxevol_gfs);

  // Step 5. Set eta damping
  if( eta_damping != -1 ) {
    griddata.params.eta_damping = eta_damping;
  }
  else {
    if( domain_size==1e6 && sinh_width==0.07 && foci_position==5 ) {
      conformally_flat_BBH_set_optimum_eta_damping(dt, CFL_FACTOR, &griddata.params);
    }
    else {
      CCTK_VERROR("Automatic computation of the optimum eta_damping parameter is only available if domain_size = 1e6 (got %e), sinh_width = 0.07 (got %e) and foci_position = 5 (got %e). Please set NRPyEllipticET::eta_damping manually on your parfile or adjust your grid parameters.",domain_size, sinh_width, foci_position);
    }
  }

  // Step 6: Print more useful information
  CCTK_INFO("Grid information:");
  CCTK_VINFO("    Points in xx0    =  %d",Nxx[0]);
  CCTK_VINFO("    Points in xx1    =  %d",Nxx[1]);
  CCTK_VINFO("    Points in xx2    =  %d",Nxx[2]);
  CCTK_VINFO("    Domain size      = %10.3e",domain_size);
  CCTK_VINFO("    sinh width       = %10.3e",sinh_width);

  CCTK_INFO("Hyperbolic relaxation information:");
  CCTK_VINFO("    Damping strength = %10.3e",griddata.params.eta_damping);
  CCTK_VINFO("    Target residual  = %10.3e",pow(10.0,log_target_residual));
  CCTK_VINFO("    Time step        = %10.3e",dt);
  CCTK_VINFO("    Final time       = %10.3e",number_of_LCT*domain_size);
  CCTK_VINFO("    Max iterations   =  %d"   ,N_final);

  // Step 7: Relaxation time loop
  // Step 7.a: Start timing for NRPyEllipticET
  const clock_t start = clock();

  // Step 7.b: Set frequency with which the L2-norm of residual is computed
  int output_l2_norm_residual_every_N = info_output_freq; //= (int)((REAL)N_final/1000.0);
  if( output_l2_norm_residual_every_N == 0 ) output_l2_norm_residual_every_N = 1;

  // Step 7.b.i: The leapfrog integrator needs uu^{-1} = uu^0 (vv = 0 initially)
  const int use_leapfrog = CCTK_EQUALS(time_integrator, "leapfrog");
  if( use_leapfrog ) {
    memcpy(griddata.gridfuncs.k_even_gfs, griddata.gridfuncs.y_n_gfs, sizeof(REAL)*NUM_EVOL_GFS*Nxx_plus_2NGHOSTS_tot);
  }
  CCTK_VINFO("    Time integrator  = %s", time_integrator);

  // Step 7.b.ii: EXPERIMENTAL grid sequencing (leapfrog only): relax on coarser grids first, then
  //              finish on this grid. See NRPyEllipticET_conformally_flat_BBH_sequencing.c.
  if( coarse_grid_levels > 0 && !use_leapfrog ) {
    CCTK_ERROR("coarse_grid_levels > 0 requires time_integrator = \"leapfrog\".");
  }
  int n = 0;
  if( coarse_grid_levels > 0 ) {
    griddata_struct coarser;
    for(int level=coarse_grid_levels; level>=1; level--) {
      const int Nc[3] = { Nxx[0]>>level, Nxx[1]>>level, Nxx[2]>>level };
      for(int d=0;d<3;d++) {
        if( Nxx[d] % (1<<(level+1)) != 0 ) {
          CCTK_VERROR("coarse_grid_levels = %d needs N0, N1, N2 divisible by %d (got %d %d %d).",
                      (int)coarse_grid_levels, 1<<(level+1), Nxx[0], Nxx[1], Nxx[2]);
        }
      }
      griddata_struct g;
      REAL dtc;
      conformally_flat_BBH_sequencing_griddata_setup(&griddata.params, Nc, CFL_FACTOR, griddata.params.eta_damping, &g, &dtc);
      const int is_coarsest = (level == coarse_grid_levels);
      if( !is_coarsest ) {
        conformally_flat_BBH_sequencing_prolongate(&coarser, &g);
        conformally_flat_BBH_sequencing_griddata_free(&coarser);
      }
      const clock_t t0 = clock();
      REAL log_res;
      const int its = conformally_flat_BBH_sequencing_relax(&g, dtc, is_coarsest ? 0 : prolongation_smoothing_steps,
                                                            coarse_log_target_residual, N_final, residual_integration_radius, &log_res);
      if( its < 0 ) CCTK_VERROR("Coarse level %dx%dx%d did not converge (log10 residual %.2f).", Nc[0], Nc[1], Nc[2], log_res);
      CCTK_VINFO("Grid sequencing: %dx%dx%d: %d iterations, log10 residual %.2f, %.2f s",
                 Nc[0], Nc[1], Nc[2], its, log_res, (clock()-t0)/((REAL)CLOCKS_PER_SEC));
      coarser = g;
    }
    conformally_flat_BBH_sequencing_prolongate(&coarser, &griddata);
    conformally_flat_BBH_sequencing_griddata_free(&coarser);
    const clock_t t0 = clock();
    REAL log_res;
    const int its = conformally_flat_BBH_sequencing_relax(&griddata, dt, prolongation_smoothing_steps, log_target_residual,
                                                          N_final, residual_integration_radius, &log_res);
    if( its < 0 ) CCTK_VERROR("Relaxation did not converge (log10 residual %.2f).", log_res);
    CCTK_VINFO("Grid sequencing: %dx%dx%d: %d iterations, %.2f s", Nxx[0], Nxx[1], Nxx[2], its, (clock()-t0)/((REAL)CLOCKS_PER_SEC));
    CCTK_VINFO("Target log(residual) (%.2lf) has been met: %.2lf after %d iterations",log_target_residual, log_res, its);
    n = 10*N_final; // skip the standard loop below
  }

  // Precomputed weights for the residual's L2 norm (bitwise identical, much cheaper per call).
  REAL *restrict l2_weights = conformally_flat_BBH_L2_norm_weights_malloc(&griddata.params, residual_integration_radius, griddata.xx);
  if( l2_weights == NULL ) CCTK_ERROR("Could not allocate L2-norm weights.");

  for(;n<=N_final;n++) { // Main loop to progress forward in time.

    // Step 7.c: Step forward one timestep (t -> t+dt) in time using
    //           chosen RK-like MoL timestepping algorithm
    //           (leapfrog also stores the residual of the previous
    //           time level in diagnostic_output_gfs as a byproduct)
    if( use_leapfrog ) {
      conformally_flat_BBH_leapfrog_step_forward_in_time(&griddata, dt);
    }
    else {
      conformally_flat_BBH_MoL_step_forward_in_time(&griddata, dt);
    }

    // Step 7.d: Check convergence & print information
    if( (n%output_l2_norm_residual_every_N == 0) || (n==N_final) ) {

      // Step 7.d.i: Compute residual and store it at diagnostic_output_gfs(UUFG)
      if( !use_leapfrog ) {
        conformally_flat_BBH_residual_all_points(&griddata.params, &griddata.rfmstruct, griddata.gridfuncs.auxevol_gfs,
                                                 griddata.gridfuncs.y_n_gfs, griddata.gridfuncs.diagnostic_output_gfs);
      }

      const int residual_gf_index   = UUGF;
      const REAL log_l2_norm_residual = conformally_flat_BBH_L2_norm_of_gf_weighted(&griddata.params, residual_gf_index, l2_weights,
                                                                                    griddata.gridfuncs.diagnostic_output_gfs);
      if( verbose && info_output_freq > 0 ) {
        // Step 7.d.ii: Compute time elapsed during the executation of this thorn
        const clock_t end = clock();
        REAL time_elapsed_s = (end-start)/((REAL)CLOCKS_PER_SEC);

        // Step 7.d.iii: Compute percentage of the run that has completed
        REAL completion = ((REAL)n)/((REAL)N_final);

        // Step 7.d.iv: Estimate how much time is left
        REAL time_left_s = time_elapsed_s*(1.0/completion - 1.0);

        // Step 7.d.v: Print information to the user
        CCTK_VINFO("It: %05d, t: %4.2lf, completion: %5.1lf%%, runtime: %3.0lf s, ETA: %3.0lf s, log(L2 norm H): %.2lf",
                   n, griddata.params.time, completion*100, time_elapsed_s, time_left_s, log_l2_norm_residual);
      }

      if( log_l2_norm_residual < log_target_residual ) {
        CCTK_VINFO("Target log(residual) (%.2lf) has been met: %.2lf after %d iterations",log_target_residual, log_l2_norm_residual, n);
        // We're done!
        n = 10*N_final;
        break;
      }
    }
  } // End main loop to progress forward in time.
  free(l2_weights);

  if( n==N_final ) {
    // This means the stopping criterion was not met. Error out.
    CCTK_ERROR("Solution did not reach the specified tolerance in the given number of iterations. Please increase max_iterations or the convergence criterion on the parfile.");
  }

  // Step 8: Initialize NRPyEllipticET_xx
  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) NRPyEllipticET_xx[0][i0] = griddata.xx[0][i0];
  for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) NRPyEllipticET_xx[1][i1] = griddata.xx[1][i1];
  for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++) NRPyEllipticET_xx[2][i2] = griddata.xx[2][i2];

  // Step 9: Initialize NRPyEllipticET_uu
#pragma omp parallel for
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0,
              0,Nxx_plus_2NGHOSTS1,
              0,Nxx_plus_2NGHOSTS2) {
    NRPyEllipticET_uu[IDX4S(UUGF,i0,i1,i2)] = griddata.gridfuncs.y_n_gfs[IDX4S(UUGF,i0,i1,i2)];
  }

  // Step 10: Print final information about the run
  const clock_t end = clock();
  CCTK_VINFO("Total execution time: %.2lf s", (end-start)/((REAL)CLOCKS_PER_SEC));

  // Step 11: Free all allocated memory
  conformally_flat_BBH_rfm_precompute_rfmstruct_freemem(&griddata.params, &griddata.rfmstruct);
  conformally_flat_BBH_freemem_bcstruct(&griddata.params, &griddata.bcstruct);
  conformally_flat_BBH_MoL_free_memory_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  conformally_flat_BBH_MoL_free_memory_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
}
