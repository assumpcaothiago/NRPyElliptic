#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include <sys/stat.h>
/*
 * // main() function:
 * Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
 * Step 1: Write test data to gridfunctions
 * Step 3: Apply curvilinear boundary conditions
 * Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
 * Step 5: Free all allocated memory
 */
int main(int argc, const char *argv[]) {

  griddata_struct griddata;
  set_Cparameters_to_default(&griddata.params);
  griddata.params.outer_bc_type = RADIATION_OUTER_BCS;

  // Step 0a:: Read command-line input, error out if nonconformant
  if(argc != 5 || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < NGHOSTS) {
    printf("Error: Expected four command-line arguments: ./NRPyElliptic_Playground Nx0 Nx1 Nx2 N_final,\n");
    printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions,\n");
    printf("and N_final is the total number of relaxation steps.\n");
    printf("Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
    exit(1);
  }
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }

  // Step 0c: Set free parameters, overwriting Cparameters defaults
  //          by hand or with command-line input, as desired.
#include "free_parameters.h"

  // Step 0d: Uniform coordinate grids are stored to *xx[3]
  // Step 0d.i: Set bcstruct
  {
    int EigenCoord;

    // Step 0d.ii: Set params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             Eigen-CoordSystem associated with the chosen CoordSystem;
    //             e.g., CoordSystem=SinhSpherical -> EigenCoord=Spherical
    EigenCoord = 1;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0e: Find ghostzone mappings; set up bcstruct
    bcstruct_set_up(&griddata.params, griddata.xx, &griddata.bcstruct);
    // Step 0e.i: Free allocated space for xx[][] array, as it was
    //            allocated in set_Nxx_dxx_invdx_params__and__xx()
    for(int i=0;i<3;i++) free(griddata.xx[i]);

    // Step 0f: Set params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //          chosen (non-Eigen) CoordSystem.
    EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
  }

  // Step 0g: Set timestep based on smallest proper distance between gridpoints and CFL factor
  REAL dt = find_timestep(&griddata.params, griddata.xx, CFL_FACTOR);
  //printf("# Timestep set to = %e\n",(double)dt);
                                    
  // Step 0h: Set total number of steps
  int N_final = atoi(argv[4]);
                                           
  int output_every_N = (int)((REAL)N_final/100.0);
  if(output_every_N == 0) output_every_N = 1;

  // Step 0j: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
  //              This is a limitation of the RK method. You are always welcome to declare & allocate
  //              additional gridfunctions by hand.
  if(NUM_AUX_GFS > NUM_EVOL_GFS) {
    printf("Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
    printf("       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
    exit(1);
  }

  // Step 0k: Declare struct for gridfunctions and allocate memory for y_n_gfs gridfunctions
  MoL_malloc_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  // Step 0l: Set up precomputed reference metric arrays
  // Step 0l.i: Allocate space for precomputed reference metric arrays.
  rfm_struct rfmstruct;
  rfm_precompute_rfmstruct_malloc(&griddata.params, &griddata.rfmstruct);

  // Step 0l.ii: Define precomputed reference metric arrays.
  rfm_precompute_rfmstruct_define(&griddata.params, griddata.xx, &griddata.rfmstruct);
  // Step 1a: Set up initial guess:
  initial_guess_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.y_n_gfs);

  // Step 1b: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  MoL_malloc_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  
  // Step 1c: Set up auxevol gridfunctions (psi_background and ADD_times_AUU)
  auxevol_gfs_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.auxevol_gfs);
  
  // Step 1d: Set up wavespeed gridfunction
  wavespeed_gf_all_points(&griddata.params, CFL_FACTOR, dt, griddata.xx, griddata.gridfuncs.auxevol_gfs);
  
  // Step 1.e: Set wavespeed at outer boundary (necessary for Sommerfeld BC)
  const REAL wavespeed_at_OB = compute_wavespeed_at_OB (&griddata.params, griddata.gridfuncs.auxevol_gfs);
  griddata.params.uu_wavespeed_at_OB = wavespeed_at_OB;
  griddata.params.vv_wavespeed_at_OB = wavespeed_at_OB;
                                                                
  // print wavespeeds for debugging
  printf("griddata.params.uu_wavespeed_at_OB = %.16e\n", griddata.params.uu_wavespeed_at_OB);
  printf("griddata.params.vv_wavespeed_at_OB = %.16e\n", griddata.params.vv_wavespeed_at_OB);
  
  // Step ?: Create output folders
  const int output_folder          = mkdir("./output_NRPyElliptic", 0777);
  const int output_folder_1D       = mkdir("./output_NRPyElliptic/1D", 0777);
  const int output_folder_2D       = mkdir("./output_NRPyElliptic/2D", 0777);
    
  // Step ?.a: Create file to store L2-norm of residual as function of relaxation time
  FILE *residual_l2_norm_file = fopen("./output_NRPyElliptic/residual_l2_norm.txt", "w"); fclose(residual_l2_norm_file);
  
  // Step ?.b: Output wavespeed gridfunction
  FILE *wavespeed_gf_1D_file = fopen("./output_NRPyElliptic/wavespeed_gf_1D.txt", "w"); 
  gridfunction_z_axis(&griddata.params, WAVESPEEDGF, griddata.xx, griddata.gridfuncs.auxevol_gfs, wavespeed_gf_1D_file);  
  fclose(wavespeed_gf_1D_file);
  
  // Step ?.c: Output the ADD_times_AUU gridfunction
  FILE *ADD_times_AUU_gf_1D_file = fopen("./output_NRPyElliptic/ADD_times_AUU_gf_1D.txt", "w");
  gridfunction_z_axis(&griddata.params, ADD_TIMES_AUUGF, griddata.xx, griddata.gridfuncs.auxevol_gfs, ADD_times_AUU_gf_1D_file);
  fclose(ADD_times_AUU_gf_1D_file);
 
  // Step ?.d: Output the psi_background gridfunction
  FILE *psi_background_gf_1D_file = fopen("./output_NRPyElliptic/psi_background_gf_1D.txt", "w");
  gridfunction_z_axis(&griddata.params, PSI_BACKGROUNDGF, griddata.xx, griddata.gridfuncs.auxevol_gfs, psi_background_gf_1D_file);
  fclose(psi_background_gf_1D_file);
  
  // Print puncture parameters for diagnostics purposes
  print_puncture_parameters(&griddata.params);
  

  for(int n=0;n<=N_final;n++) { // Main loop to progress forward in time.

    // Step 1a: Set current time to correct value & compute exact solution
    griddata.params.time = ((REAL)n)*dt;
    
    // Step 2: Diagnostics
      
    // Step 2.a: Compute residual and store it at diagnostic_output_gfs(UUFG)
    residual_all_points(&griddata.params, &griddata.rfmstruct, griddata.gridfuncs.auxevol_gfs,
                        griddata.gridfuncs.y_n_gfs, griddata.gridfuncs.diagnostic_output_gfs);
    // Step 2.b: compute and output L2-norm of residual
    const int residual_gf_index = UUGF;
    const REAL integration_radius = 100.;
    const REAL l2_norm_residual = L2_norm_of_gf(&griddata.params, residual_gf_index, integration_radius, griddata.xx, 
                                                griddata.gridfuncs.diagnostic_output_gfs);
    residual_l2_norm_file = fopen("./output_NRPyElliptic/residual_l2_norm.txt", "a");
    fprintf(residual_l2_norm_file, "%6d %10.4e %.17e\n", n, griddata.params.time, l2_norm_residual);
    fclose(residual_l2_norm_file);

    // Step 2.c: Output gridfunctions
    if (n%output_every_N==0){
      // Step 2.c.i: uu along z-axis
      char uu_1D_file_string[100];
      sprintf(uu_1D_file_string, "output_NRPyElliptic/1D/uu_1D_%08d.dat", n);
      FILE *uu_1D_file = fopen(uu_1D_file_string, "w");
      gridfunction_z_axis(&griddata.params, UUGF, griddata.xx, griddata.gridfuncs.y_n_gfs, uu_1D_file);
      fclose(uu_1D_file);
      
      // Step 2.c.ii: uu on xz-plane
      char uu_2D_file_string[100];
      sprintf(uu_2D_file_string, "output_NRPyElliptic/2D/uu_2D_%08d.dat", n);
      FILE *uu_2D_file = fopen(uu_2D_file_string, "w");
      gridfunction_xz_plane(&griddata.params, UUGF, griddata.xx, griddata.gridfuncs.y_n_gfs, uu_2D_file);
      fclose(uu_2D_file);
        
      // Step 2.c.ii: residual on xz-plane
      char residual_2D_file_string[100];
      sprintf(residual_2D_file_string, "output_NRPyElliptic/2D/residual_2D_%08d.dat", n);
      FILE *residual_2D_file = fopen(residual_2D_file_string, "w");
      gridfunction_xz_plane(&griddata.params, UUGF, griddata.xx, griddata.gridfuncs.diagnostic_output_gfs, residual_2D_file);
      fclose(residual_2D_file); 
    } // END if (n%output_every_N==0)

    // Step 3: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
    //         applying radiation outer boundary conditions.
    // Step 3.b: Step forward one timestep (t -> t+dt) in time using
    //           chosen RK-like MoL timestepping algorithm
    MoL_step_forward_in_time(&griddata, dt);

  } // End main loop to progress forward in time.
  
  // Step 4: Output uu after main loop
  FILE *uu_final_1D_file = fopen("./output_NRPyElliptic/uu_final_1D.txt", "w"); 
  gridfunction_z_axis(&griddata.params, UUGF, griddata.xx, griddata.gridfuncs.y_n_gfs, uu_final_1D_file);  
  fclose(uu_final_1D_file);

  // Step 5: Free all allocated memory  rfm_precompute_rfmstruct_freemem(&griddata.params, &griddata.rfmstruct);

  free(griddata.bcstruct.inner_bc_array); // Free inner bc data
  for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata.bcstruct.pure_outer_bc_array[ng]); // Free outer bc data
  MoL_free_memory_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  MoL_free_memory_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;}
