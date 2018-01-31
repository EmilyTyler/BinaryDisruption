#To test the outcome of a single encounter

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from encounters import integrateEncounter, impulseEncounter, calc_b_max
from impulse_test_encounter import impulseTestEncounter
from scipy.constants import au, parsec, giga, year
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Initialise variables
#Semi-major axis, m
a = 10.0**6.0 * au
print('a = ', a)
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Mass of perturbers
M_p = 1.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Number density of perturbers
n_p = rho/M_p

#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rms, a, m1, m2)

b = b_min
print('b = ', b)

#Impulse encounter
#notBound_imp, a_imp, e_imp = impulseEncounter(m1, m2, v_rms, b, a, e, M_p)
#print('a_imp = ', a_imp/au, 'au')

#Integrate encounter
#notBound_thr, a_thr, e_thr = integrateEncounter(m1, m2, v_rms, b, a, e, M_p)
#print('a_thr = ', a_thr/au, 'au')

notBound_thr, a_thr, e_thr, a_frac, e_diff = impulseTestEncounter(m1, m2, v_rms, b, a, e, M_p)
print('a_frac = ', a_frac)











