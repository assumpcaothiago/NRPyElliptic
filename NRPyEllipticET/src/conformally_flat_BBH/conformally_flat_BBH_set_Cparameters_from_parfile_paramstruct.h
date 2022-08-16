// Set params struct based on parfile
// Grid parameters
params.bScale               = foci_position;
params.AMAX                 = domain_size;
params.SINHWAA              = sinh_width;

// Parameters for puncture #1
params.puncture_0_bare_mass = conformally_flat_BBH_puncture_0_bare_mass;
params.puncture_0_x         = conformally_flat_BBH_puncture_0_pos[0];
params.puncture_0_y         = conformally_flat_BBH_puncture_0_pos[1];
params.puncture_0_z         = conformally_flat_BBH_puncture_0_pos[2];
params.puncture_0_P_x       = conformally_flat_BBH_puncture_0_P[0];
params.puncture_0_P_y       = conformally_flat_BBH_puncture_0_P[1];
params.puncture_0_P_z       = conformally_flat_BBH_puncture_0_P[2];
params.puncture_0_S_x       = conformally_flat_BBH_puncture_0_S[0];
params.puncture_0_S_y       = conformally_flat_BBH_puncture_0_S[1];
params.puncture_0_S_z       = conformally_flat_BBH_puncture_0_S[2];

// Parameters for puncture #2
params.puncture_1_bare_mass = conformally_flat_BBH_puncture_1_bare_mass;
params.puncture_1_x         = conformally_flat_BBH_puncture_1_pos[0];
params.puncture_1_y         = conformally_flat_BBH_puncture_1_pos[1];
params.puncture_1_z         = conformally_flat_BBH_puncture_1_pos[2];
params.puncture_1_P_x       = conformally_flat_BBH_puncture_1_P[0];
params.puncture_1_P_y       = conformally_flat_BBH_puncture_1_P[1];
params.puncture_1_P_z       = conformally_flat_BBH_puncture_1_P[2];
params.puncture_1_S_x       = conformally_flat_BBH_puncture_1_S[0];
params.puncture_1_S_y       = conformally_flat_BBH_puncture_1_S[1];
params.puncture_1_S_z       = conformally_flat_BBH_puncture_1_S[2];

// Lapse initial data parameter
params.initial_lapse        = initial_lapse;
params.lapse_exponent       = lapse_exponent_n;
