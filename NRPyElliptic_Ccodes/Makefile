CC     = gcc
CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -march=native -g -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99
#CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
INCLUDEDIRS = 
all: NRPyElliptic_Playground  MoL_malloc_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_free_memory_non_y_n_gfs.o MoL_step_forward_in_time.o bcstruct_set_up.o apply_bcs_outerradiation_and_inner.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o find_timestep.o xx_to_Cart.o set_Nxx_dxx_invdx_params__and__xx.o Cart_to_xx_and_nearest_i0i1i2.o Cart_to_xx_and_nearest_i0i1i2_global_grid_center.o rfm_precompute_rfmstruct_malloc.o rfm_precompute_rfmstruct_define.o rfm_precompute_rfmstruct_freemem.o set_Cparameters_to_default.o initial_guess_single_point.o initial_guess_all_points.o rhs_eval.o auxevol_gfs_all_points.o wavespeed_gf_all_points.o residual_all_points.o L2_norm_of_gf.o gridfunction_z_axis.o gridfunction_xz_plane.o print_puncture_parameters.o compute_wavespeed_at_OB.o main.o
MoL_malloc_y_n_gfs.o: MoL_malloc_y_n_gfs.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c MoL_malloc_y_n_gfs.c -o MoL_malloc_y_n_gfs.o

MoL_free_memory_y_n_gfs.o: MoL_free_memory_y_n_gfs.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c MoL_free_memory_y_n_gfs.c -o MoL_free_memory_y_n_gfs.o

MoL_malloc_non_y_n_gfs.o: MoL_malloc_non_y_n_gfs.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c MoL_malloc_non_y_n_gfs.c -o MoL_malloc_non_y_n_gfs.o

MoL_free_memory_non_y_n_gfs.o: MoL_free_memory_non_y_n_gfs.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c MoL_free_memory_non_y_n_gfs.c -o MoL_free_memory_non_y_n_gfs.o

MoL_step_forward_in_time.o: MoL_step_forward_in_time.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c MoL_step_forward_in_time.c -o MoL_step_forward_in_time.o

bcstruct_set_up.o: bcstruct_set_up.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c bcstruct_set_up.c -o bcstruct_set_up.o

apply_bcs_outerradiation_and_inner.o: apply_bcs_outerradiation_and_inner.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c apply_bcs_outerradiation_and_inner.c -o apply_bcs_outerradiation_and_inner.o

apply_bcs_inner_only.o: apply_bcs_inner_only.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c apply_bcs_inner_only.c -o apply_bcs_inner_only.o

apply_bcs_outerextrap_and_inner.o: apply_bcs_outerextrap_and_inner.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c apply_bcs_outerextrap_and_inner.c -o apply_bcs_outerextrap_and_inner.o

find_timestep.o: find_timestep.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c find_timestep.c -o find_timestep.o

xx_to_Cart.o: xx_to_Cart.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c xx_to_Cart.c -o xx_to_Cart.o

set_Nxx_dxx_invdx_params__and__xx.o: set_Nxx_dxx_invdx_params__and__xx.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c set_Nxx_dxx_invdx_params__and__xx.c -o set_Nxx_dxx_invdx_params__and__xx.o

Cart_to_xx_and_nearest_i0i1i2.o: Cart_to_xx_and_nearest_i0i1i2.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c Cart_to_xx_and_nearest_i0i1i2.c -o Cart_to_xx_and_nearest_i0i1i2.o

Cart_to_xx_and_nearest_i0i1i2_global_grid_center.o: Cart_to_xx_and_nearest_i0i1i2_global_grid_center.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c Cart_to_xx_and_nearest_i0i1i2_global_grid_center.c -o Cart_to_xx_and_nearest_i0i1i2_global_grid_center.o

rfm_precompute_rfmstruct_malloc.o: rfm_precompute_rfmstruct_malloc.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c rfm_precompute_rfmstruct_malloc.c -o rfm_precompute_rfmstruct_malloc.o

rfm_precompute_rfmstruct_define.o: rfm_precompute_rfmstruct_define.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c rfm_precompute_rfmstruct_define.c -o rfm_precompute_rfmstruct_define.o

rfm_precompute_rfmstruct_freemem.o: rfm_precompute_rfmstruct_freemem.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c rfm_precompute_rfmstruct_freemem.c -o rfm_precompute_rfmstruct_freemem.o

set_Cparameters_to_default.o: set_Cparameters_to_default.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c set_Cparameters_to_default.c -o set_Cparameters_to_default.o

initial_guess_single_point.o: initial_guess_single_point.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c initial_guess_single_point.c -o initial_guess_single_point.o

initial_guess_all_points.o: initial_guess_all_points.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c initial_guess_all_points.c -o initial_guess_all_points.o

rhs_eval.o: rhs_eval.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c rhs_eval.c -o rhs_eval.o

auxevol_gfs_all_points.o: auxevol_gfs_all_points.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c auxevol_gfs_all_points.c -o auxevol_gfs_all_points.o

wavespeed_gf_all_points.o: wavespeed_gf_all_points.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c wavespeed_gf_all_points.c -o wavespeed_gf_all_points.o

residual_all_points.o: residual_all_points.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c residual_all_points.c -o residual_all_points.o

L2_norm_of_gf.o: L2_norm_of_gf.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c L2_norm_of_gf.c -o L2_norm_of_gf.o

gridfunction_z_axis.o: gridfunction_z_axis.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c gridfunction_z_axis.c -o gridfunction_z_axis.o

gridfunction_xz_plane.o: gridfunction_xz_plane.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c gridfunction_xz_plane.c -o gridfunction_xz_plane.o

print_puncture_parameters.o: print_puncture_parameters.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c print_puncture_parameters.c -o print_puncture_parameters.o

compute_wavespeed_at_OB.o: compute_wavespeed_at_OB.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c compute_wavespeed_at_OB.c -o compute_wavespeed_at_OB.o

main.o: main.c free_parameters.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c main.c -o main.o

NRPyElliptic_Playground:   MoL_malloc_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_free_memory_non_y_n_gfs.o MoL_step_forward_in_time.o bcstruct_set_up.o apply_bcs_outerradiation_and_inner.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o find_timestep.o xx_to_Cart.o set_Nxx_dxx_invdx_params__and__xx.o Cart_to_xx_and_nearest_i0i1i2.o Cart_to_xx_and_nearest_i0i1i2_global_grid_center.o rfm_precompute_rfmstruct_malloc.o rfm_precompute_rfmstruct_define.o rfm_precompute_rfmstruct_freemem.o set_Cparameters_to_default.o initial_guess_single_point.o initial_guess_all_points.o rhs_eval.o auxevol_gfs_all_points.o wavespeed_gf_all_points.o residual_all_points.o L2_norm_of_gf.o gridfunction_z_axis.o gridfunction_xz_plane.o print_puncture_parameters.o compute_wavespeed_at_OB.o main.o
	$(CC)   MoL_malloc_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_free_memory_non_y_n_gfs.o MoL_step_forward_in_time.o bcstruct_set_up.o apply_bcs_outerradiation_and_inner.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o find_timestep.o xx_to_Cart.o set_Nxx_dxx_invdx_params__and__xx.o Cart_to_xx_and_nearest_i0i1i2.o Cart_to_xx_and_nearest_i0i1i2_global_grid_center.o rfm_precompute_rfmstruct_malloc.o rfm_precompute_rfmstruct_define.o rfm_precompute_rfmstruct_freemem.o set_Cparameters_to_default.o initial_guess_single_point.o initial_guess_all_points.o rhs_eval.o auxevol_gfs_all_points.o wavespeed_gf_all_points.o residual_all_points.o L2_norm_of_gf.o gridfunction_z_axis.o gridfunction_xz_plane.o print_puncture_parameters.o compute_wavespeed_at_OB.o main.o -o NRPyElliptic_Playground -lm -lgomp

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png NRPyElliptic_Playground
