import os
os.system("python setup.py build_ext --inplace")

import pstats, cProfile
import numpy as np
import encounters

from scipy.constants import G

#Initialise variables
#Semi-major axis, m
a = 0.1 * 3.086*10.0**16.0
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.008
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number of time steps
N_t = 5000
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0

#Number density of perturbers
n_p = rho/M_p

#Create time array
t = np.zeros(N_t)

#Keep track of a over time
A = np.zeros(N_t)
A[0] = a

#Keep track of e over time
es = np.zeros(N_t)
es[0] = e

cProfile.runctx("encounters.binning(v_rms, n_p, N_t, t, A, es, m1, m2, M_p)", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()