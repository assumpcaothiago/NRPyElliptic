# This module provides functions for setting up Curvilinear boundary conditions,
#     as documented in Tutorial-Start_to_Finish-Curvilinear_BCs_new_way.ipynb

# Authors: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#          Terrence Pierre Jacques

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC      # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict  # NRPy+: Core C code output module
from outputC import add_to_Cfunction_dict  # NRPy+: Core C code output module
from outputC import Cfunction    # NRPy+: Core C code output module
from outputC import indent_Ccode # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import finite_difference as fin  # NRPy+: Finite-difference module
import os, sys           # Standard Python modules for multiplatform OS-level functions
from UnitTesting.assert_equal import check_zero  # NRPy+: Checks whether an expression evaluates to zero.

_unused = par.Cparameters("int", __name__, "has_outer_boundary", 0)


# First register basic C data structures/macros inside NRPy_basic_defines.h
def NRPy_basic_defines_CurviBC_data_structures():
    return r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
typedef struct __ghostzone_map__ {
  short i0,i1,i2; // i0,i1,i2 stores values from -1 (used to indicate outer boundary)
                  // to Nxx_plus_2NGHOSTS*. We assume that grid extents beyond the
                  // limits of short (i.e., beyond about 32,000) are unlikely. This
                  // can be easily extended if needed, though.
} gz_map;

typedef struct __parity__ {
  int8_t parity[10]; // We store the 10 parity conditions in 10 int8_t integers,
                     // one for each condition. Note that these conditions can
                     // only take one of two values: +1 or -1, hence the use of
                     // int8_t, the smallest C data type.
} parity_condition;

typedef struct __inner_bc__ {
  gz_map inner_bc_dest_pt;
  gz_map inner_bc_src_pt;
  int8_t parity[10]; // We store the 10 parity conditions in 10 int8_t integers,
                     // one for each condition. Note that these conditions can
                     // only take one of two values: +1 or -1, hence the use of
                     // int8_t, the smallest C data type.
} inner_bc;

typedef struct __outer_bc__ {
  gz_map outer_bc_dest_pt;
  int8_t FACEi0,FACEi1,FACEi2; // FACEi* takes values of -1, 0, and +1 only,
                               // corresponding to MAXFACE, NUL, and MINFACE
                               // respectively.
                               // Thus int8_t (one byte each, the smallest C
                               // type) is sufficient.
} outer_bc;

typedef struct __bcstruct__ {
  outer_bc **outer; // Array of 1D arrays, of length
                    //   [NGHOSTS][num_ob_gz_pts[which_outer_ghostzone_point]]

  inner_bc **inner; // Array of 1D arrays, of length
                    //   [NGHOSTS][num_ib_gz_pts[which_inner_ghostzone_point]]

  // Arrays storing number of outer/inner boundary ghostzone points at each ghostzone,
  //   of length NGHOSTS:
  int     *num_ob_gz_pts;
  int     *num_ib_gz_pts;
} bc_struct;
"""


# Set unit-vector dot products (=parity) for each of the 10 parity condition types
def parity_conditions_symbolic_dot_products():
    parity = ixp.zerorank1(DIM=10)
    UnitVectors_inner = ixp.zerorank2()
    xx0_inbounds, xx1_inbounds, xx2_inbounds = sp.symbols("xx0_inbounds xx1_inbounds xx2_inbounds", real=True)
    for i in range(3):
        for j in range(3):
            UnitVectors_inner[i][j] = rfm.UnitVectors[i][j].subs(rfm.xx[0],xx0_inbounds).subs(rfm.xx[1],xx1_inbounds).subs(rfm.xx[2],xx2_inbounds)
    # Type 0: scalar
    parity[0] = sp.sympify(1)
    # Type 1: i0-direction vector or one-form
    # Type 2: i1-direction vector or one-form
    # Type 3: i2-direction vector or one-form
    for i in range(3):
        for Type in range(1, 4):
            parity[Type] += rfm.UnitVectors[Type-1][i]*UnitVectors_inner[Type-1][i]
    # Type 4: i0i0-direction rank-2 tensor
    # parity[4] = parity[1]*parity[1]
    # Type 5: i0i1-direction rank-2 tensor
    # Type 6: i0i2-direction rank-2 tensor
    # Type 7: i1i1-direction rank-2 tensor
    # Type 8: i1i2-direction rank-2 tensor
    # Type 9: i2i2-direction rank-2 tensor
    count = 4
    for i in range(3):
        for j in range(i, 3):
            parity[count] = parity[i+1]*parity[j+1]
            count = count + 1

    lhs_strings = []
    for i in range(10):
        lhs_strings.append("REAL_parity_array["+str(i)+"]")
    outstr = """
      // NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
      //      ten parity conditions, in given coordinate system.
      //      Needed for automatically determining sign of tensor across coordinate boundary.
      // Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
"""
    return outstr + outputC(parity, lhs_strings, filename="returnstring", params="preindent=4")


# For example, if the gridfunction name ends with "01", then (based on the table in the
# NRPy+ Jupyter notebook corresponding to this Python module) the set_parity_types()
# function below will set the parity_type of that gridfunction to 5. We can be assured
# this is a rather robust algorithm, because gri.register_gridfunctions() in grid.py
# will throw an error if a gridfunction's base name ends in an integer. This strict
# syntax was added with the express purpose of making it easier to set parity types
# based solely on the gridfunction name.
#
# After each parity type is found, we store the parity type of each gridfunction to
# const int8_t arrays evol_gf_parity and aux_gf_parity, appended to the end of
# NRPy_basic_defines.h.
def NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=True):
    # First add human-readable gridfunction aliases (grid.py) to NRPy_basic_defines dictionary,
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = gri.gridfunction_lists()

    # Step 3.b: set the parity conditions on all gridfunctions in gf_list,
    #       based on how many digits are at the end of their names
    def set_parity_types(list_of_gf_names):
        parity_type = []
        for name in list_of_gf_names:
            for gf in gri.glb_gridfcs_list:
                if gf.name == name:
                    parity_type__orig_len = len(parity_type)
                    if gf.DIM < 3 or gf.DIM > 4:
                        print("Error: Cannot currently specify parity conditions on gridfunctions with DIM<3 or >4.")
                        sys.exit(1)
                    if gf.rank == 0:
                        parity_type.append(0)
                    elif gf.rank == 1:
                        if gf.DIM == 3:
                            parity_type.append(int(gf.name[-1]) + 1)  # = 1 for e.g., beta^0; = 2 for e.g., beta^1, etc.
                        elif gf.DIM == 4:
                            parity_type.append(int(gf.name[-1]))  # = 0 for e.g., b4^0; = 1 for e.g., beta^1, etc.
                    elif gf.rank == 2:
                        if gf.DIM == 3:
                            # element of a list; a[-2] the
                            # second-to-last element, etc.
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            if idx0 == "0" and idx1 == "0":
                                parity_type.append(4)
                            elif (idx0 == "0" and idx1 == "1") or (idx0 == "1" and idx1 == "0"):
                                parity_type.append(5)
                            elif (idx0 == "0" and idx1 == "2") or (idx0 == "2" and idx1 == "0"):
                                parity_type.append(6)
                            elif idx0 == "1" and idx1 == "1":
                                parity_type.append(7)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(8)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(9)
                        elif gf.DIM == 4:
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            # g4DD00 = g_{tt} : parity type = 0
                            # g4DD01 = g_{tx} : parity type = 1
                            # g4DD02 = g_{ty} : parity type = 2
                            # g4DD0a = g_{ta} : parity type = a
                            if idx0 == "0":
                                parity_type.append(int(idx1))
                            elif idx1 == "0":
                                parity_type.append(int(idx0))
                            if idx0 == "1" and idx1 == "1":
                                parity_type.append(4)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(5)
                            elif (idx0 == "1" and idx1 == "3") or (idx0 == "3" and idx1 == "1"):
                                parity_type.append(6)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(7)
                            elif (idx0 == "2" and idx1 == "3") or (idx0 == "3" and idx1 == "2"):
                                parity_type.append(8)
                            elif idx0 == "3" and idx1 == "3":
                                parity_type.append(9)
                    if len(parity_type) == parity_type__orig_len:
                        print("Error: Could not figure out parity type for "+gf.gftype+" gridfunction: " + gf.name,gf.DIM,gf.name[-2],gf.name[-1],gf.rank)
                        sys.exit(1)
        if len(parity_type) != len(list_of_gf_names):
            print("Error: For some reason the length of the parity types list did not match the length of the gf list.")
            sys.exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)
    auxevol_parity_type = set_parity_types(auxevol_variables_list)

    # Output all gridfunctions to Ccodesrootdir/gridfunction_defines.h
    # ... then append to the file the parity type for each gridfunction.
    outstr = """
/* PARITY TYPES FOR ALL GRIDFUNCTIONS.
 * SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */
"""
    if len(evolved_variables_list) > 0:
        outstr += "static const int8_t evol_gf_parity[" + str(len(evolved_variables_list)) + "] = { "
        for i in range(len(evolved_variables_list) - 1):
            outstr += str(evol_parity_type[i]) + ", "
        outstr += str(evol_parity_type[len(evolved_variables_list) - 1]) + " };\n"

    if len(auxiliary_variables_list) > 0:
        outstr += "static const int8_t aux_gf_parity[" + str(len(auxiliary_variables_list)) + "] = { "
        for i in range(len(auxiliary_variables_list) - 1):
            outstr += str(aux_parity_type[i]) + ", "
        outstr += str(aux_parity_type[len(auxiliary_variables_list) - 1]) + " };\n"

    if len(auxevol_variables_list) > 0:
        outstr += "static const int8_t auxevol_gf_parity[" + str(len(auxevol_variables_list)) + "] = { "
        for i in range(len(auxevol_variables_list) - 1):
            outstr += str(auxevol_parity_type[i]) + ", "
        outstr += str(auxevol_parity_type[len(auxevol_variables_list) - 1]) + " };\n"

    if verbose == True:
        for i in range(len(evolved_variables_list)):
            print("Evolved gridfunction \"" + evolved_variables_list[i] + "\" has parity type " + str(
                evol_parity_type[i]) + ".")
        for i in range(len(auxiliary_variables_list)):
            print("Auxiliary gridfunction \"" + auxiliary_variables_list[i] + "\" has parity type " + str(
                aux_parity_type[i]) + ".")
        for i in range(len(auxevol_variables_list)):
            print("AuxEvol gridfunction \"" + auxevol_variables_list[i] + "\" has parity type " + str(
                auxevol_parity_type[i]) + ".")
    return outstr


# set_up__bc_gz_map_and_parity_condns() is performed using only the
# Eigen-Coordinate corresponding to the chosen reference_metric::CoordSystem.
#
# First we generate the C code needed for applying boundary conditions
# in generic coordinate systems, using the Eigen-Coordinate approach
# described in the Jupyter notebook documenting this Python module.
#
# EigenCoord_xx_to_Cart(): Output (x(x0,x1,x2), y(x0,x1,x2), z(x0,x1,x2))
# EigenCoord_Cart_to_xx(): Output (x0(x,y,z), x1(x,y,z), x2(x,y,z))
def EigenCoord_xx_to_Cart(i012suffix=""):
    # Step 1: Find the Eigen-Coordinate and set up the Eigen-Coordinate's reference metric:
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem",rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 2: Output C code for the Eigen-Coordinate mapping from xx->Cartesian:
    outstr = """    {
      // xx_to_Cart for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+r"""):
      REAL xx0 = xx[0][i0"""+i012suffix+"""];
      REAL xx1 = xx[1][i1"""+i012suffix+"""];
      REAL xx2 = xx[2][i2"""+i012suffix+"""];\n"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart[0]","xCart[1]","xCart[2]"],
            "returnstring", params="preindent=3")+"    }\n"

    # Step 3: Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
    rfm.reference_metric()

    # Step 4: Return EigenCoord xx_to_Cart C code
    return outstr


def EigenCoord_Cart_to_xx():
    # Step 1: Find the Eigen-Coordinate and set up the Eigen-Coordinate's reference metric:
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 2: Output the Eigen-Coordinate mapping from Cartesian->xx:
    # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        print("ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for")
        print("       reference_metric::CoordSystem = "+par.parval_from_str("reference_metric::CoordSystem"))
        print("       Boundary conditions in curvilinear coordinates REQUiRE this be set.")
        sys.exit(1)
    # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    outstr = """    // Cart_to_xx for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+");\n"
    outstr += outputC([rfm.Cart_to_xx[0],rfm.Cart_to_xx[1],rfm.Cart_to_xx[2]],
                      ["Cart_to_xx0_inbounds","Cart_to_xx1_inbounds","Cart_to_xx2_inbounds"],
                      filename="returnstring", params="preindent=2")

    # Step 3: Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
    rfm.reference_metric()

    # Step 4: Return EigenCoord Cart_to_xx C code
    return outstr


# Next we set up set_up__bc_gz_map_and_parity_condns(), which loops over all points
# on the numerical grid (interior and ghost zone points included), with the aim of
# filling gz_map *bc_gz_map and parity_condition *bc_parity_conditions at each point.
# To do so, the function implements the algorithm described above in Step 1.
#
# That is to say, at each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
#
# Step 1. Convert the curvilinear coordinate (x0,x1,x2) to the corresponding Cartesian
#         coordinate (x,y,z)
# Step 2. Find the Cartesian grid point (i0_inbounds,i1_inbounds,i2_inbounds) in the
#         grid interior or outer boundary corresponding to this Cartesian coordinate
#         (x,y,z).
# Step 3. If and only if we are on an outer boundary ghost zone or in the grid interior,
#         i0_inbounds==i0, i1_inbounds==i1, and i2_inbounds==i2, and inner boundary
#         conditions do not apply: set bc_gz_map to  (-1,-1,-1) , and for all 10
#         gridfunction parities, set parity=1.
# Step 4. If i0_inbounds==i0, i1_inbounds==i1, and i2_inbounds==i2, does not hold true,
#         then (i0,i1,i2) is an inner boundary point: set bc_gz_map to
#         (i0_inbounds,i1_inbounds,i2_inbounds), and for all 10 gridfunction parities,
#         evaluate all dot products of the unit vectors evaluated at
#         (i0_inbounds,i1_inbounds,i2_inbounds) and (i0,i1,i2), as prescribed above in
#         Step 1. The C code for computing the needed symbolic dot products is generated
#         above in Step 2.
def add_to_Cfunction_dict_set_up__bc_gz_map_and_parity_condns(rel_path_to_Cparams=os.path.join(".")):
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    desc = ""
    c_type = "void"
    name = "set_up__bc_gz_map_and_parity_condns"
    params = """const paramstruct *restrict params,
                                             REAL *restrict xx[3], gz_map *restrict bc_gz_map,parity_condition *restrict bc_parity_conditions"""
    body = r"""
  // xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx0;
  // -> xxmin[0] = xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0
  const REAL xxmin[3] = { xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
                          xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
                          xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };
  //fprintf(stderr,"hey inside setbc: %e %e %e | %e %e\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],dxx0);
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0,0,Nxx_plus_2NGHOSTS1,0,Nxx_plus_2NGHOSTS2) {
    // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
    REAL xCart[3];
""" + EigenCoord_xx_to_Cart(i012suffix="") + r"""
    //EigenCoord_xx_to_Cart(params, xx, i0,i1,i2, xCart);
    REAL Cartx = xCart[0];
    REAL Carty = xCart[1];
    REAL Cartz = xCart[2];

    // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
    //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
    //      the point is an outer boundary point.
    //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
    //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
    //      appropriate parity condition (+/- 1).
    REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
""" + EigenCoord_Cart_to_xx() + r"""
    const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
    const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
    const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );

    // Step 2.a: (Sanity/validation check) Convert the interior point
    //           x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) to Cartesian coordinates,
    //           make sure that the Cartesian coordinate matches the Cartesian coordinate of
    //           x0(i0),x1(i1),x2(i2). If not, error out!
    REAL xCart_orig[3]; for(int ii=0;ii<3;ii++) xCart_orig[ii] = xCart[ii];
    //EigenCoord_xx_to_Cart(params, xx, i0_inbounds,i1_inbounds,i2_inbounds, xCart);
""" + EigenCoord_xx_to_Cart(i012suffix="_inbounds") + r"""

//fprintf(stderr,"Cartesian agreement: ( %.15e %.15e %.15e ) ?= ( %.15e %.15e %.15e )\n",
// (double)xCart_orig[0],(double)xCart_orig[1],(double)xCart_orig[2],
// (double)xCart[0],(double)xCart[1],(double)xCart[2]);

#define EPS_ABS 1e-8
    if(fabs( (double)(xCart_orig[0] - xCart[0]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[1] - xCart[1]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[2] - xCart[2]) ) > EPS_ABS) {
       fprintf(stderr,"Error. """+CoordSystem+r""": Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
               (double)xCart_orig[0],(double)xCart_orig[1],(double)xCart_orig[2],
               (double)xCart[0],(double)xCart[1],(double)xCart[2],
               xx[0][i0],xx[1][i1],xx[2][i2],
               xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
               Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
      exit(1);
    }

    // Step 3: Set bc_gz_map and bc_parity_conditions.
    if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0) {
      // Step 3.a: Iff we are on an outer boundary point or in the grid
      //           interior, i0_inbounds==i0, i1_inbounds==i1, and
      //           i2_inbounds==i2, and inner boundary conditions do not
      //           apply: set bc_gz_map to -1, and parity=1.
      bc_gz_map[IDX3S(i0,i1,i2)].i0=-1;
      bc_gz_map[IDX3S(i0,i1,i2)].i1=-1;
      bc_gz_map[IDX3S(i0,i1,i2)].i2=-1;
      for(int which_parity=0; which_parity<10; which_parity++) {
        bc_parity_conditions[IDX3S(i0,i1,i2)].parity[which_parity] = 1;
      }
    } else {
      // Step 3.b: If we are on an *inner* boundary point:
      // 1. Set bc_gz_map at (i0,i1,i2) to the point
      //    in the interior to which this boundary
      //    point maps, and
      // 2. Perform the unit vector dot products
      //    necessary to set all 10 possible parity
      //    conditions, calling function
      //    set_parity_from_unit_vector_dot_product()
      bc_gz_map[IDX3S(i0,i1,i2)].i0=i0_inbounds;
      bc_gz_map[IDX3S(i0,i1,i2)].i1=i1_inbounds;
      bc_gz_map[IDX3S(i0,i1,i2)].i2=i2_inbounds;
      const REAL xx0 = xx[0][i0];
      const REAL xx1 = xx[1][i1];
      const REAL xx2 = xx[2][i2];
      const REAL xx0_inbounds = xx[0][i0_inbounds];
      const REAL xx1_inbounds = xx[1][i1_inbounds];
      const REAL xx2_inbounds = xx[2][i2_inbounds];
      REAL REAL_parity_array[10];
      {
      // Evaluate dot products needed for setting parity
      //     conditions at a given point (xx0,xx1,xx2),
      //     using C code generated by NRPy+
"""+parity_conditions_symbolic_dot_products()+r"""
      }
      //eval_symbolic_dot_products_to_set_parity_conditions(params,  REAL_parity_array,  xx0,xx1,xx2,
      //                                                          xx0_inbounds,xx1_inbounds,xx2_inbounds);
      for(int whichparity=0;whichparity<10;whichparity++) {
          //printf("Good? Parity %d evaluated to %e\n",whichparity,(double)REAL_parity_array[whichparity]);
          // Perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
          if( (REAL_parity_array[whichparity]  > 0 && fabs(REAL_parity_array[whichparity] - (+1)) > 1e-8) ||
              (REAL_parity_array[whichparity] <= 0 && fabs(REAL_parity_array[whichparity] - (-1)) > 1e-8) ) {
              fprintf(stderr,"Error at point (%d %d %d); (%e %e %e); maps to (%e %e %e).\n",
                      i0,i1,i2, xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
              fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
                      REAL_parity_array[whichparity]);
              exit(1);
          }
          if(REAL_parity_array[whichparity] < 0.0) bc_parity_conditions[IDX3S(i0,i1,i2)].parity[whichparity] = -1;
          if(REAL_parity_array[whichparity] > 0.0) bc_parity_conditions[IDX3S(i0,i1,i2)].parity[whichparity] = +1;
      } // END for(int whichparity=0;whichparity<10;whichparity++)
    } // END if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0)
  } // END LOOP_REGION(0,Nxx_plus_2NGHOSTS0,0,Nxx_plus_2NGHOSTS1,0,Nxx_plus_2NGHOSTS2)
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)


# As described above, set_up__bc_gz_map_and_parity_condns() sets gz_map *bc_gz_map
# and parity_condition *bc_parity_conditions at all grid points. While this
# information could be used directly to apply boundary conditions,
# it is more efficient (both in memory and CPU time) to instead store only the
# needed information to the bcstruct array, so that when applying boundary
# conditions, we can simply loop over bcstruct.
#
# The algorithm follows the bc_struct data type defined above in Step 2,
# and loops over all boundary ghost zones, from the innermost layer (which_gz=0)
# outward (to which_gz=NGHOSTS-1). This is necessary, as, for example, some
# ghost zone points on the which_gz=1 layer depend on ghost zones being set on the
# which_gz=0 layer.
def add_to_Cfunction_dict_set_bcstruct():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]

    prefunc = """
static const int8_t MAXFACE = -1;
static const int8_t NUL     = +0;
static const int8_t MINFACE = +1;
"""
    desc = """
set_bcstruct() loops from the innermost boundary
     ghostzones on the cube ("which_gz==0",
     corresponding to the single layer of ghostzones
     closest to the interior data), and at each
     ghostzone layer, we apply the following 5-step
     algorithm:
Step 1: Count the number of outer and inner
        boundary points, store to
        num_ob_pts and num_ib_pts, respectively.
Step 2: Now that we know the number of outer
        boundary points on this ghostzone layer,
        allocate memory needed for storing the
        outer and inner boundary condition data.
Step 2.a: At all outer boundary ghost zones, allocate
          memory for a single member of the outer_bc
          data type.
Step 2.b: At all inner boundary ghost zones, allocate
          memory for a single member of the inner_bc
          data type.
Step 3: Store the number of outer and inner boundary
        points on each ghostzone layer, where e.g.,
        which_gz==0 corresponds to the innermost
        ghostzones on the numerical domain.
Step 4: Store information needed for outer boundary
        conditions, to outer_bc_dest_pt and
        outer_bc_face arrays.
Step 5: Store information needed for inner boundary
        conditions, including interior point to which
        inner ghost zone maps, and parity conditions
        for all 10 gridfunction types.
"""
    c_type = "void"

    name = "set_bcstruct"
    params = """const paramstruct *restrict params,
                  gz_map *restrict bc_gz_map,
                  parity_condition *bc_parity_conditions,
                  bc_struct *restrict bcstruct"""
    body = r"""
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };

  // Loop from the innermost ghostzone on the cube (which_gz==0) and work outward.
  //      This ordering is necessary, as ghostzones at which_gz==1 will generally
  //      depend on ghostzones at which_gz==0 being already set.
  for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {

    // Step 1: Count the number of outer and inner
    //         boundary points, store to
    //         num_ob_pts and num_ib_pts, respectively.
"""
    body += r"""
#define COUNT_INNER_OR_OUTER if(bc_gz_map[IDX3S(i0,i1,i2)].i0==-1) { num_ob_pts++;} else { num_ib_pts++; }
"""
    body += r"""
    int num_ob_pts = 0;
    int num_ib_pts = 0;
    LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imin[0]--;
    LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imax[0]++;
    LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imin[1]--;
    LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imax[1]++;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) { COUNT_INNER_OR_OUTER } imin[2]--;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) { COUNT_INNER_OR_OUTER } imax[2]++;
    
    // Step 2: Now that we know the number of outer boundary points on this ghostzone
    //    layer, we allocate memory needed for storing the outer and inner boundary
    //     condition data.
    
    // Step 2.a: At all outer boundary ghost zones, allocate memory for a single member of the outer_bc
    //           data type.
    bcstruct->outer[which_gz] = (outer_bc *)malloc(sizeof(outer_bc)*num_ob_pts);
    // Step 2.b: At all inner boundary ghost zones, allocate memory for a single member of the inner_bc
    //           data type.
    bcstruct->inner[which_gz] = (inner_bc *)malloc(sizeof(inner_bc)*num_ib_pts);
    
    // Step 3: Store the number of outer and inner boundary points on each ghostzone layer, where e.g.,
    //         which_gz==0 corresponds to the innermost ghostzones on the numerical domain.
    bcstruct->num_ob_gz_pts[which_gz] = num_ob_pts;
    bcstruct->num_ib_gz_pts[which_gz] = num_ib_pts;
    
    // Reset imin[] and imax[], to prepare for the next step.
    for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;}
    
    // Step 4: Store information needed for outer boundary conditions, to outer_bc_dest_pt[which_gz][]
//         and outer_bc_face[which_gz][] arrays:
"""
    body += """#define OB_SET(facei0,facei1,facei2) if(bc_gz_map[IDX3S(i0,i1,i2)].i0==-1) { \\"""

    body += r"""
          bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0 = i0;       \
          bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1 = i1;       \
          bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2 = i2;       \
          bcstruct->outer[which_gz][pt].FACEi0= facei0;                 \
          bcstruct->outer[which_gz][pt].FACEi1= facei1;                 \
          bcstruct->outer[which_gz][pt].FACEi2= facei2;                 \
          pt++; }

    int pt = 0;
    LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) {OB_SET(MINFACE,NUL,NUL)} imin[0]--;
    LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) {OB_SET(MAXFACE,NUL,NUL)} imax[0]++;
    LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) {OB_SET(NUL,MINFACE,NUL)} imin[1]--;
    LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) {OB_SET(NUL,MAXFACE,NUL)} imax[1]++;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) {OB_SET(NUL,NUL,MINFACE)} imin[2]--;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) {OB_SET(NUL,NUL,MAXFACE)} imax[2]++;
    // fprintf(stderr,"num OB points with which_gz = %d: %d | should be: %d\n",which_gz,pt,num_ob_gz_pts[which_gz]);
    
    // Reset imin[] and imax[], to prepare for the next step.
    for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;}
    
    // Step 5: Store information needed for inner boundary conditions, including interior point to which
    //         inner ghost zone maps, and parity conditions for all 10 gridfunction types.
#define IB_SET if(bc_gz_map[IDX3S(i0,i1,i2)].i0!=-1) {                  \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0=i0;               \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1=i1;               \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2=i2;               \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0 =bc_gz_map[IDX3S(i0,i1,i2)].i0; \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1 =bc_gz_map[IDX3S(i0,i1,i2)].i1; \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2 =bc_gz_map[IDX3S(i0,i1,i2)].i2; \
    for(int ii=0;ii<10;ii++) {                                          \
    bcstruct->inner[which_gz][pt].parity[ii] =                          \
      (int8_t)bc_parity_conditions[IDX3S(i0,i1,i2)].parity[ii]; }       \
    pt++; }

    pt = 0;
    LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) {IB_SET} imin[0]--;
    LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) {IB_SET} imax[0]++;
    LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) {IB_SET} imin[1]--;
    LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) {IB_SET} imax[1]++;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) {IB_SET} imin[2]--;
    LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) {IB_SET} imax[2]++;

  } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
"""
    add_to_Cfunction_dict(
        includes=includes, prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body),
        rel_path_to_Cparams=os.path.join("."))


# driver_bcstruct() allocates memory for boundary condition pointers
# and calls functions referenced above, including
# set_up__bc_gz_map_and_parity_condns() and set_bcstruct(). It then
# frees unneeded memory after bcstruct has been fully initialized.
def add_to_Cfunction_dict_driver_bcstruct(enable_mask=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """driver_bcstruct(): Set up bcstruct.
WARNING: Needs Nxx_plus_2NGHOSTS_tot * [sizeof(gz_map) + sizeof(parity_condition) ] bytes of temporary memory!
   (This memory will be deallocated at end of function, in freemem_bcstruct().)
   Note that gz_map consists of 3*sizeof(short = 2 bytes), and parity_condition consists of 10*sizeof(int8_t = 1 byte)
   Thus the total cost is about 16 bytes at all gridpoints, or 2 double-precision gridfunctions.
   To avoid memory overload, be sure to set this up prior to allocation of other gridfunctions.
STEPS:
1. Allocate memory for bc_gz_map, which maps inner gzs to the
   appropriate interior gridpoint. In the case of outer gzs, the
   gz point maps to itself.
2. Allocate storage for parity_condition, which
"""
    c_type = "void"
    name = "driver_bcstruct"
    params = """const paramstruct *restrict params, bc_struct *restrict bcstruct, REAL *restrict xx[3]"""
    set_bcstruct_extra_arg = ""
    body = """
const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

// Step 1: Allocate memory storage for bc_gz_map, which
//         in the case a boundary point is a *parity*
//         boundary, is set to the interior, non-
//         boundary point corresponding to the same
//         Cartesian gridpoint. Otherwise bc_gz_map
//         is set to (i0,i1,i2) = (-1,-1,-1).
gz_map *restrict bc_gz_map = (gz_map *restrict)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);

// Step 2: Allocate memory storage for bc_parity_conditions,
//         which store parity conditions for all 10
//         gridfunction types at all grid points.
parity_condition *restrict bc_parity_conditions = (parity_condition *restrict)malloc(sizeof(parity_condition)*Nxx_plus_2NGHOSTS_tot);

// Step 3: Set bc_gz_map and bc_parity_conditions at *all*
//         points; on the boundary and otherwise.
set_up__bc_gz_map_and_parity_condns(params, xx, bc_gz_map,
                                    bc_parity_conditions);

// Step 4: Declare and allocate memory for bcstruct,
//         which will store all information needed for
//         applying the boundary conditions.
bcstruct->outer = (outer_bc **)malloc(sizeof(outer_bc *)*NGHOSTS);
bcstruct->inner = (inner_bc **)malloc(sizeof(inner_bc *)*NGHOSTS);
bcstruct->num_ob_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);
bcstruct->num_ib_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);

// Step 5: Store all information needed to quickly and
//         efficiently apply boundary conditions. This
//         function transfers all information from
//         bc_gz_map (defined at *all gridpoints*) into
//         bcstruct (defined only at boundary points).
//         Thus when this function has finished,
//         bc_gz_map is no longer needed.
set_bcstruct(params,bc_gz_map,
             bc_parity_conditions,
             bcstruct"""+set_bcstruct_extra_arg+r""");
// Do not apply outer boundary conditions on grids that do not possess the outer boundary!
//if(params->has_outer_boundary == 0) {
//  for(int i=0;i<NGHOSTS;i++) bcstruct->num_ob_gz_pts[i] = 0;
//}

// Step 6: As described in Step 4, bc_gz_map is no
//         longer needed at this point, so we free its
//         memory. Farewell, friend!
free(bc_gz_map);
free(bc_parity_conditions);
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=os.path.join("."))


# Free memory allocated within bcstruct.
def add_to_Cfunction_dict_freemem_bcstruct():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = "Free memory allocated within bcstruct"
    c_type = "void"
    name = "freemem_bcstruct"
    params = "const paramstruct *restrict params, const bc_struct *restrict bcstruct"
    body = r"""  for(int i=0;i<NGHOSTS;i++) { free(bcstruct->outer[i]);  free(bcstruct->inner[i]); }
  free(bcstruct->outer);  free(bcstruct->inner);
  free(bcstruct->num_ob_gz_pts); free(bcstruct->num_ib_gz_pts);
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=os.path.join("."))


###############################
## RADIATION (NewRad-like) BOUNDARY CONDITION FUNCTIONS
def get_arb_offset_FD_coeffs_indices(FDORDER, offset, deriv):
    # deriv = 1 <-- 1st derivative
    Minv = fin.setup_FD_matrix__return_inverse(FDORDER+1, offset)
    indices = []
    coeffs = []
    for i in range(FDORDER+1):
        indices.append(i-int(FDORDER/2) + offset)
        coeffs.append(Minv[i, deriv])
    return coeffs, indices

def setup_Cfunction_FD1_arbitrary_upwind(dirn, BC_FDORDER=-1):
    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if BC_FDORDER == -1:
        BC_FDORDER = default_FDORDER

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", BC_FDORDER)

    includes = []
    desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
    c_type = "static inline REAL"
    name = "FD1_arbitrary_upwind_x"+str(dirn)+"_dirn"
    params = """const paramstruct *restrict params, const REAL *restrict gf,
    const int i0,const int i1,const int i2, const int offset"""
    body = r"""switch(offset) {
"""
    tmp_list = []
    for offset in range(0, int(BC_FDORDER / 2) + 1):
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += "case " + str(offset) + ":\n"
        body += "  return ("
        coeffs, indices = get_arb_offset_FD_coeffs_indices(BC_FDORDER, offset, 1)
        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue  # skip this iteration if coeff=0
            offset = str(indices[i])
            if i > 0:
                body += "          "
            if offset == "0":
                body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0+"+offset+",i1,i2)]\n"
                elif dirn == 1:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1+"+offset+",i2)]\n"
                elif dirn == 2:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2+"+offset+")]\n"
        body = body[:-1].replace("+-", "-") + ") * invdx"+str(dirn)+";\n"
    body += """}
return 0.0 / 0.0;  // poison output if offset computed incorrectly
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", default_FDORDER)

    return func


def compute_Jacobian_and_inverseJacobian_tofrom_Spherical():
    # Step 2.a: First construct Jacobian matrix:
    Jac_dUSph_dDrfmUD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            Jac_dUSph_dDrfmUD[i][j] = sp.diff(rfm.xxSph[i], rfm.xx[j])
    Jac_dUrfm_dDSphUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSph_dDrfmUD)
    return Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD


def setup_Cfunction_compute_partial_r_f(BC_FDORDER=-1):
    desc = "Compute \partial_r f"
    c_type = "static inline REAL"
    name = "compute_partial_r_f"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs,
                                       const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
                                       const int FACEi0,const int FACEi1,const int FACEi2,
                                       const REAL partial_x0_r, const REAL partial_x1_r, const REAL partial_x2_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()

    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if BC_FDORDER == -1:
        BC_FDORDER = default_FDORDER

    body = r"""  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = BC_FDORDER/2 = """ + str(int(BC_FDORDER/2)) + r"""
  const int FD1_stencil_radius = """ + str(int(BC_FDORDER/2)) + r""";

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {i0,i1,i2}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
"""
    for i in range(3):
        si = str(i)
        if check_zero(Jac_dUrfm_dDSphUD[i][0]):
            body += "  const REAL partial_x"+si+"_f=0.0;\n"
        else:
            body += "  int i"+si+"_offset = FACEi"+si+";  // up/downwind on the faces. This offset should never go out of bounds.\n"
            body += "  if(dest_i"+si+" < FD1_stencil_radius) i"+si+"_offset = FD1_stencil_radius-dest_i"+si+";\n"
            body += "  else if(dest_i"+si+" > (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1)) i"+si+"_offset = (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1) - dest_i"+si+";\n"
            body += "  const REAL partial_x"+si+"_f=FD1_arbitrary_upwind_x"+si+"_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i"+si+"_offset);\n\n"
    body += r"""  return partial_x0_r*partial_x0_f + partial_x1_r*partial_x1_f + partial_x2_r*partial_x2_f;
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func


def setup_Cfunction_r_and_partial_xi_r_derivs():
    desc = "Compute r(xx0,xx1,xx2)."
    c_type = "static inline void"
    name = "r_and_partial_xi_r_derivs"
    params = """const paramstruct *restrict params,const REAL xx0,const REAL xx1,const REAL xx2,
                                  REAL *r, REAL *rinv, REAL *partial_x0_r,REAL *partial_x1_r,REAL *partial_x2_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    body = outputC([rfm.xxSph[0], 1/rfm.xxSph[0],
                    Jac_dUrfm_dDSphUD[0][0], Jac_dUrfm_dDSphUD[1][0], Jac_dUrfm_dDSphUD[2][0]],
                   ["*r", "*rinv", "*partial_x0_r", "*partial_x1_r", "*partial_x2_r"], filename="returnstring",
                   params="preindent=1,outCverbose=False,includebraces=False")
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func


def setup_Cfunction_radiation_bcs_curvilinear(BC_FDORDER=-1):
    includes = []
    prefunc = ""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    for i in range(3):
        # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
        if not check_zero(Jac_dUrfm_dDSphUD[i][0]):
            prefunc += setup_Cfunction_FD1_arbitrary_upwind(dirn=i, BC_FDORDER=BC_FDORDER)
    prefunc += setup_Cfunction_r_and_partial_xi_r_derivs()
    prefunc += setup_Cfunction_compute_partial_r_f(BC_FDORDER=BC_FDORDER)
    desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
    c_type = "static inline void"
    name = "radiation_bcs_curvilinear"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct,REAL *restrict xx[3],
                           const REAL *restrict gfs, REAL *restrict gfs_rhss,
                           const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
                           const int FACEi0,const int FACEi1,const int FACEi2"""
    body = r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
REAL r,rinv, partial_x0_r,partial_x1_r,partial_x2_r;
REAL r_int,r_intinv, partial_x0_r_int,partial_x1_r_int,partial_x2_r_int;
r_and_partial_xi_r_derivs(params,xx[0][dest_i0],    xx[1][dest_i1],    xx[2][dest_i2],    &r,    &rinv,    &partial_x0_r,    &partial_x1_r,    &partial_x2_r);
r_and_partial_xi_r_derivs(params,xx[0][dest_i0_int],xx[1][dest_i1_int],xx[2][dest_i2_int],&r_int,&r_intinv,&partial_x0_r_int,&partial_x1_r_int,&partial_x2_r_int);
const REAL partial_r_f     = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_r,    partial_x1_r,    partial_x2_r);
const REAL partial_r_f_int = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_r_int,partial_x1_r_int,partial_x2_r_int);

const int idx3 = IDX3S(dest_i0,dest_i1,dest_i2);
const int idx3_int = IDX3S(dest_i0_int,dest_i1_int,dest_i2_int);

const REAL partial_t_f_int = gfs_rhss[IDX4ptS(which_gf, idx3_int)];

const REAL c = gridfunctions_wavespeed[which_gf];
const REAL f_infinity = gridfunctions_f_infinity[which_gf];
const REAL f = gfs[IDX4ptS(which_gf, idx3)];
const REAL f_int = gfs[IDX4ptS(which_gf, idx3_int)];
const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) * r_intinv);

const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

gfs_rhss[IDX4ptS(which_gf, idx3)] = partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=includes, prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func
###############################


# apply_bcs_curvilinear() loops over all NUM_GFS gridfunctions in the gfs
# IDX4 array and, using a bcstruct filled in by the set_bcstruct() function
# above, applies boundary conditions to each ghost zone layer, starting
# with the innermost layer and working outward.
#
# This function is meant to be called within a Method of Lines timestepping
# algorithm, at or near the end of each substep.
def add_to_Cfunction_dict_apply_bcs_curvilinear(outer_bcs_type="extrapolation", inner_only=False, BC_FDORDER=-1):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    # Set prefunc:
    if not inner_only and outer_bcs_type == "extrapolation":
        prefunc = r"""// Declare boundary condition EXTRAP_BC_UPDATE_OUTER macro,
//          which updates a single outer boundary point
//          of the 3D grid cube using quadratic polynomial
//          extrapolation.
#define EXTRAP_BC_UPDATE_OUTER(which_gf, i0,i1,i2, FACEX0,FACEX1,FACEX2) { \
    const int idx3 = IDX3S(i0,i1,i2);                                   \
    gfs[IDX4S(which_gf,i0,i1,i2)] =                                     \
      +3.0*gfs[IDX4S(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]     \
      -3.0*gfs[IDX4S(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]     \
      +1.0*gfs[IDX4S(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)];    \
  }
"""
    elif not inner_only and outer_bcs_type == "radiation":
        prefunc = setup_Cfunction_radiation_bcs_curvilinear(BC_FDORDER=BC_FDORDER)
    elif not inner_only:
        print("outer_bcs_type == " + outer_bcs_type + " NOT SUPPORTED.")
        sys.exit(1)
    else:  # if inner_only:
        prefunc = ""
    desc = r"""Curvilinear boundary condition driver routine: Apply BCs to all six
  boundary faces of the 3D numerical domain, filling in the
  innermost ghost zone layer first, and moving outward.
"""
    c_type = "void"
    name = "apply_bcs_curvilinear_" + outer_bcs_type
    if inner_only:
        name = "apply_bcs_curvilinear_inner_only"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict xx[3],
                           REAL *restrict gfs"""
    if not inner_only and outer_bcs_type == "radiation":
        params += ", REAL *restrict gfs_rhss"
    body = r"""#pragma omp parallel for
  for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
"""
    if not inner_only:
        body += r"""
      // First apply OUTER boundary conditions,
      //   in case an INNER (parity) boundary point
      //   needs data at the outer boundary:
      // After updating each face, adjust imin[] and imax[]
      //   to reflect the newly-updated face extents.
      for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
"""
        if outer_bcs_type == "radiation":
            body += r"""        // *** Apply radiation BCs to all outer boundary points. ***
        radiation_bcs_curvilinear(params, bcstruct, xx, gfs, gfs_rhss,  which_gf,
                                  bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0,
                                  bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1,
                                  bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2,
                                  bcstruct->outer[which_gz][pt].FACEi0,
                                  bcstruct->outer[which_gz][pt].FACEi1,
                                  bcstruct->outer[which_gz][pt].FACEi2);
"""
        elif outer_bcs_type == "extrapolation":
            body += r"""        // *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
        EXTRAP_BC_UPDATE_OUTER(which_gf,
                               bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0,
                               bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1,
                               bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2,
                               bcstruct->outer[which_gz][pt].FACEi0,
                               bcstruct->outer[which_gz][pt].FACEi1,
                               bcstruct->outer[which_gz][pt].FACEi2);
"""
        body += """      }\n\n"""
    body += r"""      // Apply INNER (parity) boundary conditions:
      for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
        const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
        const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
        const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
        const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
        const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
        const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
"""
    inner_bc_str = """        gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
          bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src,i1src,i2src)];"""
    if outer_bcs_type == "radiation":
        body += inner_bc_str.replace("gfs[IDX", "gfs_rhss[IDX")
    else:
        body += inner_bc_str
    body += r"""
      } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
    } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
  } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
"""
    add_to_Cfunction_dict(
        includes=includes, prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=os.path.join("."))


# Only call this after ALL gridfunctions have been registered!
def CurviBoundaryConditions_register_C_functions_and_NRPy_basic_defines(verbose=True):
    # First register C functions needed by CurviBCs
    add_to_Cfunction_dict_set_bcstruct()
    add_to_Cfunction_dict_driver_bcstruct()
    # DEPEND ON COORDSYSTEM:
    # add_to_Cfunction_dict_set_up__bc_gz_map_and_parity_condns()
    # add_to_Cfunction_dict_apply_bcs_curvilinear()
    add_to_Cfunction_dict_apply_bcs_curvilinear(inner_only=True)

    add_to_Cfunction_dict_freemem_bcstruct()
    # Then set up the dictionary entry for CurviBC in NRPy_basic_defines
    Nbd_str  = NRPy_basic_defines_CurviBC_data_structures()
    Nbd_str += NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=verbose)
    outC_NRPy_basic_defines_h_dict["CurviBoundaryConditions"] = Nbd_str
