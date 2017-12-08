#To test the outcome of a single encounter

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from encounters import encounter
from encounters import calc_b_max

#Initialise variables
#Semi-major axis, m
a = 0.1 * 3.086*10.0**16.0
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Mass of perturbers
M_p = 1000.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p

#Minimum impact parameter
b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rms, a, m1, m2)

'''
b = b_max
v = v_rms
(notBound, a_new, e_new) = encounter(m1, m2, v, b, a, e, M_p)

print('Checking b_max')
print('a_frac = ', (a_new-a)/a)
'''

a = 18.0 * 1.5*10.0**11.0
b = b_min
v = v_rms
(notBound, a_new, e_new) = encounter(m1, m2, v, b, a, e, M_p)

print('Checking binary breaking')
print('Binary broken?', notBound)
