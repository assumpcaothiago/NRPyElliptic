// Set params struct based on parfile
// Grid parameters
griddata.params.bScale               = foci_position;
griddata.params.AMAX                 = domain_size;
griddata.params.SINHWAA              = sinh_width;

// Parameters for puncture #1
griddata.params.puncture_0_bare_mass = conformally_flat_BBH_puncture_0_bare_mass;
griddata.params.puncture_0_x         = conformally_flat_BBH_puncture_0_pos[0];
griddata.params.puncture_0_y         = conformally_flat_BBH_puncture_0_pos[1];
griddata.params.puncture_0_z         = conformally_flat_BBH_puncture_0_pos[2];
griddata.params.puncture_0_P_x       = conformally_flat_BBH_puncture_0_P[0];
griddata.params.puncture_0_P_y       = conformally_flat_BBH_puncture_0_P[1];
griddata.params.puncture_0_P_z       = conformally_flat_BBH_puncture_0_P[2];
griddata.params.puncture_0_S_x       = conformally_flat_BBH_puncture_0_S[0];
griddata.params.puncture_0_S_y       = conformally_flat_BBH_puncture_0_S[1];
griddata.params.puncture_0_S_z       = conformally_flat_BBH_puncture_0_S[2];

// Parameters for puncture #2
griddata.params.puncture_1_bare_mass = conformally_flat_BBH_puncture_1_bare_mass;
griddata.params.puncture_1_x         = conformally_flat_BBH_puncture_1_pos[0];
griddata.params.puncture_1_y         = conformally_flat_BBH_puncture_1_pos[1];
griddata.params.puncture_1_z         = conformally_flat_BBH_puncture_1_pos[2];
griddata.params.puncture_1_P_x       = conformally_flat_BBH_puncture_1_P[0];
griddata.params.puncture_1_P_y       = conformally_flat_BBH_puncture_1_P[1];
griddata.params.puncture_1_P_z       = conformally_flat_BBH_puncture_1_P[2];
griddata.params.puncture_1_S_x       = conformally_flat_BBH_puncture_1_S[0];
griddata.params.puncture_1_S_y       = conformally_flat_BBH_puncture_1_S[1];
griddata.params.puncture_1_S_z       = conformally_flat_BBH_puncture_1_S[2];

// Lapse initial data parameter
griddata.params.initial_lapse        = initial_lapse;
griddata.params.lapse_exponent       = lapse_exponent_n;
