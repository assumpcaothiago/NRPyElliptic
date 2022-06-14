
// Free parameters related to physical system:
griddata.params.time = 0.0; // Initial simulation time corresponds to time=0.

// Free parameters related to numerical timestep:
REAL CFL_FACTOR = 0.7;

// Set free-parameter values.

const REAL domain_size    = 1000000.0;
const REAL sinh_width     = 0.07;
const REAL sinhv2_const_dr= 0.05;
const REAL SymTP_bScale   = 5.0;

griddata.params.bScale =  SymTP_bScale;
griddata.params.AMAX   =  domain_size;
        griddata.params.SINHWAA = sinh_width;


// Free parameters related to punctures
griddata.params.eta_damping = 11.0;
griddata.params.puncture_0_bare_mass = 0.5184199353358704;
griddata.params.puncture_1_bare_mass = 0.39193567996522616;
griddata.params.puncture_0_x = 0.0;
griddata.params.puncture_0_y = 0.0;
griddata.params.puncture_0_z = 5.0;
griddata.params.puncture_1_x = 0.0;
griddata.params.puncture_1_y = 0.0;
griddata.params.puncture_1_z = -5.0;
griddata.params.puncture_0_x = 0.09530152296974252;
griddata.params.puncture_0_y = 0.0;
griddata.params.puncture_0_z = -0.00084541526517121;
griddata.params.puncture_1_P_x = -0.09530152296974252;
griddata.params.puncture_1_P_y = 0.0;
griddata.params.puncture_1_P_z = 0.00084541526517121;
griddata.params.puncture_0_S_x = 0.0;
griddata.params.puncture_0_S_y = 0.09509112426035504;
griddata.params.puncture_0_S_z = 0.0;
griddata.params.puncture_1_S_x = 0.0;
griddata.params.puncture_1_S_y = -0.09156449704142013;
griddata.params.puncture_1_S_z = 0.0;
