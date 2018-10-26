#To test the impulse approximation

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import au, parsec, giga, year

from impulse_test_encounter import impulseTestEncounter, impulseTestEquations
from encounters import calc_b_max
from random_direction import randomDirection
from random_binary import setupRandomBinary

from internal_units import *

print('l_scale =', length_scale())
print('m_scale =', mass_scale())
print('t_scale =', time_scale())

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0 / mass_scale()
m2 = 2.0*10.0**30.0 / mass_scale()
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0 / mass_scale()
#RMS of Maxwellian velocity distribution, m/s
v_rel = 1.0 * 220.0 * 1000.0 * time_scale()/length_scale()
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
rho = rho /mass_scale() * length_scale()**(3.0)
#Number density of perturbers
n_p = rho/M_p

#Semi-major axis
a = 10.0**5.0 * au / length_scale()
#Minimum impact parameter
b_min = (np.pi*n_p*v_rel*(10.0*giga*year/time_scale()))**(-0.5)
print('b_min=', b_min)
#Maximum impact parameter for impulse approximation
b_max = v_rel * np.sqrt(a**3.0/(G()*(m1+m2)))
print('b_max=', b_max)
#Impact parameter
b = 10.0**4.0 * au /length_scale()
print('b =', b)

#Number of encounters
N_enc = 10

#
E_frac_avg = 0.0

for i in range(N_enc):
        notBound, a_thr, e_thr, a_frac, e_diff, E_frac = impulseTestEncounter(m1, m2, v_rel, b, a, e, M_p)
        E_frac_avg += E_frac
        
#Normalise
E_frac_avg /= N_enc

print('Average fractional error in energy change =', E_frac_avg)
