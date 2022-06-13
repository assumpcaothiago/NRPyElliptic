const int Nxx0 = params.Nxx0;  // grid::Nxx0
const int Nxx1 = params.Nxx1;  // grid::Nxx1
const int Nxx2 = params.Nxx2;  // grid::Nxx2
const int Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0;  // grid::Nxx_plus_2NGHOSTS0
const int Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1;  // grid::Nxx_plus_2NGHOSTS1
const int Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2;  // grid::Nxx_plus_2NGHOSTS2
const REAL dxx0 = params.dxx0;  // grid::dxx0
const REAL dxx1 = params.dxx1;  // grid::dxx1
const REAL dxx2 = params.dxx2;  // grid::dxx2
const REAL xxmin0 = params.xxmin0;  // grid::xxmin0
const REAL xxmin1 = params.xxmin1;  // grid::xxmin1
const REAL xxmin2 = params.xxmin2;  // grid::xxmin2
const REAL xxmax0 = params.xxmax0;  // grid::xxmax0
const REAL xxmax1 = params.xxmax1;  // grid::xxmax1
const REAL xxmax2 = params.xxmax2;  // grid::xxmax2
const REAL invdx0 = params.invdx0;  // grid::invdx0
const REAL invdx1 = params.invdx1;  // grid::invdx1
const REAL invdx2 = params.invdx2;  // grid::invdx2
const REAL Cart_originx = params.Cart_originx;  // grid::Cart_originx
const REAL Cart_originy = params.Cart_originy;  // grid::Cart_originy
const REAL Cart_originz = params.Cart_originz;  // grid::Cart_originz
const REAL Cart_CoM_offsetx = params.Cart_CoM_offsetx;  // grid::Cart_CoM_offsetx
const REAL Cart_CoM_offsety = params.Cart_CoM_offsety;  // grid::Cart_CoM_offsety
const REAL Cart_CoM_offsetz = params.Cart_CoM_offsetz;  // grid::Cart_CoM_offsetz
const REAL bScale = params.bScale;  // reference_metric::bScale
const REAL SINHWAA = params.SINHWAA;  // reference_metric::SINHWAA
const REAL AMAX = params.AMAX;  // reference_metric::AMAX
const REAL eta_damping = params.eta_damping;  // NRPyElliptic_codegen.NRPyElliptic_RHSs::eta_damping
const REAL time = params.time;  // NRPyElliptic_codegen.NRPyElliptic_RHSs::time
const REAL bare_mass_0 = params.bare_mass_0;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::bare_mass_0
const REAL bare_mass_1 = params.bare_mass_1;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::bare_mass_1
const REAL puncture_0_x = params.puncture_0_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_0_x
const REAL puncture_0_y = params.puncture_0_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_0_y
const REAL puncture_0_z = params.puncture_0_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_0_z
const REAL puncture_1_x = params.puncture_1_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_1_x
const REAL puncture_1_y = params.puncture_1_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_1_y
const REAL puncture_1_z = params.puncture_1_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::puncture_1_z
const REAL P0_x = params.P0_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P0_x
const REAL P0_y = params.P0_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P0_y
const REAL P0_z = params.P0_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P0_z
const REAL P1_x = params.P1_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P1_x
const REAL P1_y = params.P1_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P1_y
const REAL P1_z = params.P1_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::P1_z
const REAL S0_x = params.S0_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S0_x
const REAL S0_y = params.S0_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S0_y
const REAL S0_z = params.S0_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S0_z
const REAL S1_x = params.S1_x;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S1_x
const REAL S1_y = params.S1_y;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S1_y
const REAL S1_z = params.S1_z;  // NRPyElliptic_codegen.NRPyElliptic_SourceTerm::S1_z
const int has_outer_boundary = params.has_outer_boundary;  // CurviBoundaryConditions.CurviBoundaryConditions_new_way::has_outer_boundary
