#Evolves a population of binaries
import numpy as np
import os
os.system("python setup.py build_ext --inplace")
from monte_carlo import MCEncounters

#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0**10.0*365.25*24.0*60.0*60.0

#Initial binary population
N_bin = 10

(a_final, e_final) = MCEncounters(v_rms, n_p, T, m1, m2, M_p, a_initial, e_initial, N_bin)

