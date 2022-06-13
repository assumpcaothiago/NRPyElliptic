# Common parameters for Puncture ID evolutions
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Thiago Assumpcao
#          assumpcaothiago **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par    # NRPy+: Parameter interface

thismodule = __name__

# Parameters common to/needed by all ScalarWave Python modules

# Step P2a: Define the C parameter wavespeed. The `wavespeed`
#           variable is a proper SymPy variable, so it can be
#           used in below expressions. In the C code, it acts
#           just like a usual parameter, whose value is
#           specified in the parameter file.
wavespeed = par.Cparameters("REAL", thismodule, "wavespeed",1.0)
# Step P2a: Define the C parameter eta_damping.
eta_damping = par.Cparameters("REAL", thismodule, "eta_damping",1.0)

# Step P3: Define the C parameters for punctures.
# Step P3a: Define bare masses
bare_mass_0 = par.Cparameters("REAL",thismodule,"bare_mass_0", 1.0)
bare_mass_1 = par.Cparameters("REAL",thismodule,"bare_mass_1", 0.0)
# Step P3b.i: Define puncture 0 positions
puncture_0_x = par.Cparameters("REAL",thismodule,"puncture_0_x", 0.0)
puncture_0_y = par.Cparameters("REAL",thismodule,"puncture_0_y", 0.0)
puncture_0_z = par.Cparameters("REAL",thismodule,"puncture_0_z", 0.0)
puncture_0 = [puncture_0_x, puncture_0_y, puncture_0_z] # define list with punctures
# Step P3b.ii: Define puncture 1 positions
puncture_1_x = par.Cparameters("REAL",thismodule,"puncture_1_x", 0.0)
puncture_1_y = par.Cparameters("REAL",thismodule,"puncture_1_y", 0.0)
puncture_1_z = par.Cparameters("REAL",thismodule,"puncture_1_z", 0.0)
puncture_1 = [puncture_1_x, puncture_1_y, puncture_1_z] # define list with punctures
# Step P3c.i: Define linear momentum 0
P0_x = par.Cparameters("REAL",thismodule,"P0_x", 0.0)
P0_y = par.Cparameters("REAL",thismodule,"P0_y", 0.0)
P0_z = par.Cparameters("REAL",thismodule,"P0_z", 0.0)
P0 = [P0_x, P0_y, P0_z] # define list with linear momenta
# Step P3c.ii: Define linear momentum 1
P1_x = par.Cparameters("REAL",thismodule,"P1_x", 0.0)
P1_y = par.Cparameters("REAL",thismodule,"P1_y", 0.0)
P1_z = par.Cparameters("REAL",thismodule,"P1_z", 0.0)
P1 = [P1_x, P1_y, P1_z] # define list with linear momenta
# Step P3d.i: Define angular momentum 0
S0_x = par.Cparameters("REAL",thismodule,"S0_x", 0.0)
S0_y = par.Cparameters("REAL",thismodule,"S0_y", 0.0)
S0_z = par.Cparameters("REAL",thismodule,"S0_z", 0.5)
S0 = [S0_x, S0_y, S0_z] # define list with angular momenta
# Step P3d.ii: Define angular momentum 1
S1_x = par.Cparameters("REAL",thismodule,"S1_x", 0.0)
S1_y = par.Cparameters("REAL",thismodule,"S1_y", 0.0)
S1_z = par.Cparameters("REAL",thismodule,"S1_z", 0.0)
S1 = [S1_x, S1_y, S1_z] # define list with angular momenta
