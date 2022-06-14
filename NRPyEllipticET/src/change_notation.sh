#!/bin/bash

# Puncture 0
sed -i 's/bare_mass_0/puncture_0_bare_mass/g' *.c *.h
sed -i 's/P0_x/puncture_0_P_x/g' *.c *.h
sed -i 's/P0_y/puncture_0_P_y/g' *.c *.h
sed -i 's/P0_z/puncture_0_P_z/g' *.c *.h
sed -i 's/S0_x/puncture_0_S_x/g' *.c *.h
sed -i 's/S0_y/puncture_0_S_y/g' *.c *.h
sed -i 's/S0_z/puncture_0_S_z/g' *.c *.h

# Puncture 1
sed -i 's/bare_mass_1/puncture_1_bare_mass/g' *.c *.h
sed -i 's/P1_x/puncture_1_P_x/g' *.c *.h
sed -i 's/P1_y/puncture_1_P_y/g' *.c *.h
sed -i 's/P1_z/puncture_1_P_z/g' *.c *.h
sed -i 's/S1_x/puncture_1_S_x/g' *.c *.h
sed -i 's/S1_y/puncture_1_S_y/g' *.c *.h
sed -i 's/S1_z/puncture_1_S_z/g' *.c *.h
