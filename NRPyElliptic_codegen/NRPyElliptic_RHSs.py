# Generating C code for the right-hand-side
#  of the scalar wave equation, in
#  ***curvilinear*** coordinates, using a
#  reference-metric formalism
#
# Authors: Thiago Assumpcao
#          assumpcaothiago **at** gmail **dot* com
#
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-PunctureInitialDataCurvilinear.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import grid as gri                            # NRPy+: Functionality for handling numerical grids
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm                # NRPy+: Reference metric support
#import PunctureInitialData.CommonParams       # NRPy+: Common parameters for all Puncture Initial Data modules
import sympy as sp

# The name of this module ("PunctureInitialDataCurvilinear") is given by __name__:
thismodule = __name__


def PunctureInitialDataCurvilinear_RHSs():
    # Step 1: Get the spatial dimension, defined in the
    #         NRPy+ "grid" module. With reference metrics,
    #         this must be set to 3 or fewer.
    DIM = par.parval_from_str("DIM")

    # Step 2: Set up the reference metric and
    #         quantities derived from the
    #         reference metric.
    rfm.reference_metric()
    
    ############################################
    # This only works with azimuthally-symmetric data
    #par.set_parval_from_str("indexedexp::symmetry_axes",'2')
    ############################################
    
    # Step 3: Define the C parameters. The variables are proper
    #         SymPy variables, so thet can be used in below expressions.
    #         In the C code, they act just like a usual parameter,
    #         whose value is specified in the parameter file.         
    
    # Wave speed is a gridfunction, so the parameter wavespeed is not defined    
    #wavespeed = par.Cparameters("REAL",thismodule,"wavespeed", 1.0) 
    eta_damping = par.Cparameters("REAL",thismodule,"eta_damping", 1.0)
    time = par.Cparameters("REAL", thismodule, "time", 0.0)
    # Declare wavespeeds of uu and vv at outer boundary
    uu_wavespeed_at_OB = par.Cparameters("REAL", thismodule, "uu_wavespeed_at_OB", 1.0)
    vv_wavespeed_at_OB = par.Cparameters("REAL", thismodule, "vv_wavespeed_at_OB", 1.0)
    
    # Declare spatially-dependent wavespeed as grid function
    wavespeed = gri.register_gridfunctions("AUXEVOL", ["wavespeed"])

    # Step 4: Define psi_background and ADD_times_AUU
    psi_background, ADD_times_AUU = gri.register_gridfunctions("AUXEVOL",["psi_background","ADD_times_AUU"])

    # Step 5: Compute the contracted Christoffel symbols:
    contractedGammahatU = ixp.zerorank1()
    for k in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                contractedGammahatU[k] += rfm.ghatUU[i][j] * rfm.GammahatUDD[k][i][j]

    # Step 6: Register gridfunctions that are needed as input
    #         to the scalar wave RHS expressions.
    uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])

    # Step 7a: Declare the rank-1 indexed expression \partial_{i} u,
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    uu_dD = ixp.declarerank1("uu_dD")

    # Step 7b: Declare the rank-2 indexed expression \partial_{ij} u,
    #          which is symmetric about interchange of indices i and j
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    uu_dDD = ixp.declarerank2("uu_dDD","sym01")

    # Step 8: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)
    global uu_rhs,vv_rhs

    # Step 9: Define right-hand sides for the evolution.
    # Step 9a: uu_rhs = v - eta_damping*uu:
    uu_rhs = vv - eta_damping*uu
    # Step 9b: The right-hand side of the \partial_t v equation
    #          is given by:
    #          \hat{g}^{ij} \partial_i \partial_j u - \hat{\Gamma}^i \partial_i u.
    #          ^^^^^^^^^^^^ PART 1 ^^^^^^^^^^^^^^^^ ^^^^^^^^^^ PART 2 ^^^^^^^^^^^
    vv_rhs = sp.sympify(0)
    for i in range(DIM):
        # PART 2:
        vv_rhs -= contractedGammahatU[i]*uu_dD[i]
        for j in range(DIM):
            # PART 1:
            vv_rhs += rfm.ghatUU[i][j]*uu_dDD[i][j]
    
    # Add source term
    vv_rhs += sp.Rational(1,8)*ADD_times_AUU*sp.Pow(psi_background + uu, -7)
    
    # Before multiplying vv_rhs by wavespeed**2, set residual
    global residual
    residual = vv_rhs
    
    # Update RHS 
    vv_rhs *= wavespeed*wavespeed


    # Step 10: Generate C code for scalarwave evolution equations,
    #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])

 
                                                        
