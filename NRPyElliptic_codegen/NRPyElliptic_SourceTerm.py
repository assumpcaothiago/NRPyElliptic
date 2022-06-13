# Generating C code for the source term
#  of the Hamiltonian constraint equation, in
#  ***curvilinear*** coordinates, using a
#  reference-metric formalism
#
# Author: Thiago Assumpcao
#         assumpcaothiago **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-PunctureInitialDataCurvilinear_RHSs.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import grid as gri                            # NRPy+: Functionality for handling numerical grids
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm                # NRPy+: Reference metric support
#import PunctureInitialData.CommonParams       # NRPy+: Common parameters for all Puncture Initial Data modules
import sympy as sp

# The name of this module ("PunctureInitialDataCurvilinear") is given by __name__:
thismodule = __name__

# Step 2.a: Define dot and cross product of vectors
def dot(vec1,vec2):
    vec1_dot_vec2 = sp.sympify(0)
    for i in range(3):
        vec1_dot_vec2 += vec1[i]*vec2[i]
    return vec1_dot_vec2
def cross(vec1,vec2):
    vec1_cross_vec2 = ixp.zerorank1()
    LeviCivitaSymbol = ixp.LeviCivitaSymbol_dim3_rank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                vec1_cross_vec2[i] += LeviCivitaSymbol[i][j][k]*vec1[j]*vec2[k]
    return vec1_cross_vec2

# define function to replace xxCart by the Cartesian coordinates expressed in curvilinear coordinates
def replace_cart_coord_by_xx (src, DIM=3):
    for k in range(DIM):
        src = src.subs(rfm.Cart[k],rfm.xx_to_Cart[k])
    return src

# Step (?): Define function to compute vector field V^i
def VU_cart_two_punctures(x0U, P0U, S0U, x1U, P1U, S1U):
    """
    Evaluate VU with numerical values or symbolic expressions

    xnU = position of n-th puncture
    PnU = linear momentum of n-th puncture
    SnU = spin of n-th puncture

    """
    DIM = par.parval_from_str("DIM")
    
    def VU_cart_single_puncture(xU, PU, SU):
        r = sp.sympify(0)
        for i in range(DIM):
            r += sp.Pow(xU[i], 2)
        r = sp.sqrt(r)

        VU_cart = ixp.zerorank1()
        for i in range(DIM):
            VU_cart[i] += sp.Rational(-7, 4) * PU[i] / r \
            + sp.Rational(-1, 4) * dot(xU, PU) * xU[i] / sp.Pow(r, 3) \
            + cross(xU, SU)[i] / sp.Pow(r, 3)  
        return VU_cart

    V0U_cart = VU_cart_single_puncture(x0U, P0U, S0U)
    V1U_cart = VU_cart_single_puncture(x1U, P1U, S1U)

    VU_cart = ixp.zerorank1()
    for i in range(DIM):   
        VU_cart[i] = V0U_cart[i] + V1U_cart[i]     
    return VU_cart

# Step (?): Define function to compute the components of the tracefree conformal
#           extrinsic curvature in Cartesian coordinates with two lower indices
def ADD_conf_cartesian(VU_cart):

    """
    Input: vector VU in Cartesian coordinates

    Output: Components of the conformally-related extrinsic 
            curvature in Cartesian coordinates with two lower indices
    """
    DIM = par.parval_from_str("DIM")
    
    # Step (?): Define Cartesian grid
    xxCart = rfm.Cart

    # Step 1: Define vector field V_i
    #           In Cartesian coordinates, V_i = V^i
    VD_cart = ixp.zerorank1()
    for i in range(DIM):
        VD_cart[i] = VU_cart[i]

    # Step 2: Compute the derivatives: sp.diff(VD_cart[i], xxCart[j])
    VD_cart_dD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            VD_cart_dD[i][j] = sp.diff(VD_cart[i], xxCart[j])

    # Step 3: Compute the divergence
    V_cart_div = sp.sympify(0)
    for i in range(DIM):
        V_cart_div += VD_cart_dD[i][i]

    # Step 4: Define Kronecker delta symbol as a function
    def kronecker_delta(i,j):
        if (i==j):
            return sp.sympify(1)
        else:
            return sp.sympify(0)

    # Step 5: Compute the components of the conformally-related extrinsic curvature
    #           in Cartesian coordinates with two lower indices, and return result
    ADD_conf_cart = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            ADD_conf_cart[i][j] = (VD_cart_dD[i][j] + VD_cart_dD[j][i]
                                   - sp.Rational(2,3)*kronecker_delta(i,j)*V_cart_div)
    return ADD_conf_cart

# Step (?): Define function to compute the contraction of the tracefree conformal
#           extrinsic curvature in Cartesian coordinates
def ADD_times_AUU_conf_cartesian(ADD_conf_cart):
    """
    Input: Components of the conformally-related extrinsic 
            curvature in Cartesian coordinates with two lower indices

    Output: The contraction ADD_conf_cart[i][j]*ADD_conf_cart[i][j]
    """
    # Set grid dimension
    DIM = par.parval_from_str("DIM")
    ADD_times_AUU_conf_cart = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            ADD_times_AUU_conf_cart += ADD_conf_cart[i][j]*ADD_conf_cart[i][j]
    return ADD_times_AUU_conf_cart

# Define function to return numeric values of ADD_times_AUU
def ADD_times_AUU_conf_cartesian_numeric(punc0U, P0U, S0U, punc1U, P1U, S1U, xyz_array):
    """
    Inputs:
        - puncnU    = positions of n-th puncture (3-D vector)
        - PnU       = linear momentum of of n-th puncture (3-D vector)
        - SnU       = angular momentum of n-th puncture (3-D vector)
        - xyz_array = array of 3-D vectors with Cartesian components of each point
                      where ADD_times_AUU is to be calculated
    Output:
        - ADD_times_AUU = array of length = len(xyz_array) with the numerical value
                          of ADD_times_AUU at each point in xyz_array
    """
    # Set grid dimension
    DIM = par.parval_from_str("DIM")
    # Declare relative position of punctures symbolically
    xU0 = [ (rfm.Cart[i]-punc0U[i]) for i in range(DIM)]
    xU1 = [ (rfm.Cart[i]-punc1U[i]) for i in range(DIM)]
    # Compute ADD_times_AUU symbolically
    VU_cart = VU_cart_two_punctures(xU0, P0U, S0U, xU1, P1U, S1U)
    ADD_conf_cart = ADD_conf_cartesian(VU_cart)
    ADD_times_AUU_conf_cart = ADD_times_AUU_conf_cartesian(ADD_conf_cart)
    # Create list where to store ADD_times_AUU at each point in xyz_array
    ADD_times_AUU_list = []
    # Replace symbolic value of Cartesian coordinates by xyz
    for xyz in xyz_array:
        ADD_times_AUU = ADD_times_AUU_conf_cart
        for i in range(3):
            ADD_times_AUU = ADD_times_AUU.subs(rfm.Cart[i], xyz[i])
        ADD_times_AUU_list.append(ADD_times_AUU)
    return ADD_times_AUU_list


# Step (?): Define function to compute the background conformal 
#           factor in Cartesian coordinates
def psi_background_cartesian (bare_mass_0, xU0, bare_mass_1, xU1):
    """
    Inputs: 
        - bare mass of each puncture
        - xnU = relative position of n-th puncture

    Output: background conformal factor in Cartesian coordinates
    """
    # Set grid dimension
    DIM = par.parval_from_str("DIM")
    
    def psi_background_cartesian_single_puncture (bare_mass, xU):
        # Compute distance of puncture from the origin
        r = sp.sympify(0)
        for i in range(DIM):
            r += sp.Pow(xU[i], 2)
        r = sp.sqrt(r)
        
        return (bare_mass/(2*r))
    
    psi_puncture_0 = psi_background_cartesian_single_puncture (bare_mass_0, xU0)
    psi_puncture_1 = psi_background_cartesian_single_puncture (bare_mass_1, xU1)
    
    return (1 + psi_puncture_0 + psi_puncture_1)

def psi_background_cartesian_numeric (bare_mass_0, punc0U, bare_mass_1, punc1U, xyz_array):
    """
    Inputs:
        - bare_mass_n = bare mass of the n-th puncture
        - puncnU      = positions of n-th puncture (3-D vector)
        - xyz_array = array of 3-D vectors with Cartesian components of each point
                      where psi_background is to be calculated
    Output:
        - psi_background = array of length = len(xyz_array) with the numerical value
                          of psi_background at each point in xyz_array
    """
    # Set grid dimension
    DIM = par.parval_from_str("DIM")
    # Declare relative position of punctures symbolically
    xU0 = [ (rfm.Cart[i]-punc0U[i]) for i in range(DIM)]
    xU1 = [ (rfm.Cart[i]-punc1U[i]) for i in range(DIM)]
    # Compute psi_background symbolically
    psi_background_sym = psi_background_cartesian (bare_mass_0, xU0, bare_mass_1, xU1)
    # Create list where to store psi_background at each point in xyz_array
    psi_background_list = []
    # Replace symbolic value of Cartesian coordinates by xyz
    for xyz in xyz_array:
        psi_background = psi_background_sym
        for i in range(DIM):
            psi_background = psi_background.subs(rfm.Cart[i], xyz[i])
        psi_background_list.append(psi_background)
    return psi_background_list


def PunctureInitialDataCurvilinear_psi_background_and_ADD_AUU_term():
    # Step 0a: Define the C parameters. The variables are 
    #           proper SymPy variables, so they can be
    #           used in below expressions. In the C code, they act
    #           just like usual parameters, whose value are
    #           specified in the parameter file.
    #wavespeed = par.Cparameters("REAL", thismodule, "wavespeed",1.0)
    # Step 0a: Define the C parameter eta_damping.
    #eta_damping = par.Cparameters("REAL", thismodule, "eta_damping",1.0)
    # Step 0: Define the C parameters for punctures.
    # Step 0a: Define bare masses
    bare_mass_0 = par.Cparameters("REAL",thismodule,"bare_mass_0", 1.0)
    bare_mass_1 = par.Cparameters("REAL",thismodule,"bare_mass_1", 0.0)
    # Step 0b.i: Define puncture 0 positions
    puncture_0_x = par.Cparameters("REAL",thismodule,"puncture_0_x", 0.0)
    puncture_0_y = par.Cparameters("REAL",thismodule,"puncture_0_y", 0.0)
    puncture_0_z = par.Cparameters("REAL",thismodule,"puncture_0_z", 0.0)
    puncture_0 = [puncture_0_x, puncture_0_y, puncture_0_z] # define list with punctures
    # Step 0b.ii: Define puncture 1 positions
    puncture_1_x = par.Cparameters("REAL",thismodule,"puncture_1_x", 0.0)
    puncture_1_y = par.Cparameters("REAL",thismodule,"puncture_1_y", 0.0)
    puncture_1_z = par.Cparameters("REAL",thismodule,"puncture_1_z", 0.0)
    puncture_1 = [puncture_1_x, puncture_1_y, puncture_1_z] # define list with punctures
    # Step 0c.i: Define linear momentum 0
    P0_x = par.Cparameters("REAL",thismodule,"P0_x", 0.0)
    P0_y = par.Cparameters("REAL",thismodule,"P0_y", 0.0)
    P0_z = par.Cparameters("REAL",thismodule,"P0_z", 0.0)
    P0 = [P0_x, P0_y, P0_z] # define list with linear momenta
    # Step 0c.ii: Define linear momentum 1
    P1_x = par.Cparameters("REAL",thismodule,"P1_x", 0.0)
    P1_y = par.Cparameters("REAL",thismodule,"P1_y", 0.0)
    P1_z = par.Cparameters("REAL",thismodule,"P1_z", 0.0)
    P1 = [P1_x, P1_y, P1_z] # define list with linear momenta
    # Step 0d.i: Define angular momentum 0
    S0_x = par.Cparameters("REAL",thismodule,"S0_x", 0.0)
    S0_y = par.Cparameters("REAL",thismodule,"S0_y", 0.0)
    S0_z = par.Cparameters("REAL",thismodule,"S0_z", 0.2)
    S0 = [S0_x, S0_y, S0_z] # define list with angular momenta
    # Step 0d.ii: Define angular momentum 1
    S1_x = par.Cparameters("REAL",thismodule,"S1_x", 0.0)
    S1_y = par.Cparameters("REAL",thismodule,"S1_y", 0.0)
    S1_z = par.Cparameters("REAL",thismodule,"S1_z", 0.0)
    S1 = [S1_x, S1_y, S1_z] # define list with angular momenta

    # Step 1: Get the spatial dimension, defined in the
    #         NRPy+ "grid" module. With reference metrics,
    #         this must be set to 3 or fewer.
    DIM = par.parval_from_str("DIM")

    # Step 2: Set up the reference metric and
    #         quantities derived from the
    #         reference metric.
    rfm.reference_metric()
    
    # Step (?): Define grid
    xx = gri.xx
    # Step (?): Define Cartesian grid
    xxCart = rfm.Cart
#     # Step (?): Register the Cartesian grid as a NRPy+ AUX gridfunction
#     Cartx, Carty, Cartz = gri.register_gridfunctions("AUX",["Cartx", "Carty", "Cartz"])
                
    # Step (?): Define linear momenta vectors
    PU = ixp.zerorank2()
    for i in range(DIM):
        PU[0][i] = P0[i]
        PU[1][i] = P1[i]
        
    # Step (?): Define angular momenta vectors
    SU = ixp.zerorank2()
    for i in range(DIM):
        SU[0][i] = S0[i]
        SU[1][i] = S1[i]
        
    # Step (?): Set relative positions of each puncture
    # xU = ixp.declarerank2("xU", "nosym")
    xU = ixp.zerorank2()
    for i in range(DIM):
        xU[0][i] = xxCart[i] - puncture_0[i]
        xU[1][i] = xxCart[i] - puncture_1[i]

    # Step (?): Define vector field V_i
    VU_cart = VU_cart_two_punctures(xU[0], PU[0], SU[0], xU[1], PU[1], SU[1])
        
    # Step (?): Define extrinsic curvature
    ADD_conf_cart = ADD_conf_cartesian(VU_cart)
    
    # Step (?): Compute the scalar ADD*AUU in Cartesian coordinates
    ADD_times_AUU_conf_cart = ADD_times_AUU_conf_cartesian(ADD_conf_cart)
    
    # Step(?): Compute psi_background
    psi_background_cart = psi_background_cartesian (bare_mass_0, xU[0], bare_mass_1, xU[1])
    global psi_background
    psi_background = replace_cart_coord_by_xx(psi_background_cart)
    
    # Step(?): Compute ADD_times_AUU
    global ADD_times_AUU
    ADD_times_AUU = replace_cart_coord_by_xx(ADD_times_AUU_conf_cart)