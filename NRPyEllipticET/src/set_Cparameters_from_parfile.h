// Set params struct based on parfile
// Grid parameters
griddata.params.bScale               = foci_position;
griddata.params.AMAX                 = domain_size;
griddata.params.SINHWAA              = sinh_width;

// Parameters for puncture #1
griddata.params.puncture_0_bare_mass = puncture_0_bare_mass;
griddata.params.puncture_0_x         = puncture_0_pos[0];
griddata.params.puncture_0_y         = puncture_0_pos[1];
griddata.params.puncture_0_z         = puncture_0_pos[2];
griddata.params.puncture_0_P_x       = puncture_0_P[0];
griddata.params.puncture_0_P_y       = puncture_0_P[1];
griddata.params.puncture_0_P_z       = puncture_0_P[2];
griddata.params.puncture_0_S_x       = puncture_0_S[0];
griddata.params.puncture_0_S_y       = puncture_0_S[1];
griddata.params.puncture_0_S_z       = puncture_0_S[2];

// Parameters for puncture #2
griddata.params.puncture_1_bare_mass = puncture_1_bare_mass;
griddata.params.puncture_1_x         = puncture_1_pos[0];
griddata.params.puncture_1_y         = puncture_1_pos[1];
griddata.params.puncture_1_z         = puncture_1_pos[2];
griddata.params.puncture_1_P_x       = puncture_1_P[0];
griddata.params.puncture_1_P_y       = puncture_1_P[1];
griddata.params.puncture_1_P_z       = puncture_1_P[2];
griddata.params.puncture_1_S_x       = puncture_1_S[0];
griddata.params.puncture_1_S_y       = puncture_1_S[1];
griddata.params.puncture_1_S_z       = puncture_1_S[2];

// Lapse initial data parameter
griddata.params.initial_lapse        = initial_lapse;
griddata.params.lapse_exponent       = lapse_exponent_n;

__attribute__ ((unused)) REAL CFL_FACTOR = 0.7;
