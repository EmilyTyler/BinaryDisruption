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
v_rel = 220.0 * 1000.0 * time_scale()/length_scale()
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
print('b_min, au =', b_min*length_scale()/au)
#Maximum impact parameter for impulse approximation
b_max = v_rel * np.sqrt(a**3.0/(G()*(m1+m2)))
print('b_max, au =', b_max*length_scale()/au)

#Impact parameter
b = 10.0**2.5 * au / length_scale()
print('b, au =', b*length_scale()/au)
#Number of encounters
N_enc = 10**5
print('N_enc =', N_enc)
#Simulation parameters
delta=10.0**(-4.0)
eta=0.02

#
E_frac_avg = 0.0
E_frac_var = 0.0
E_frac_error_avg = 0.0
E_frac_error_var = 0.0

E_ini = np.zeros(N_enc)
E_thr = np.zeros(N_enc)
E_imp = np.zeros(N_enc)
bs = np.zeros(N_enc)

for i in range(N_enc):
        notBound, a_thr, e_thr, a_frac, e_diff, E_frac, E_frac_error, E_ini[i], E_thr[i], E_imp[i], bs[i] = impulseTestEncounter(m1, m2, v_rel, b, a, e, M_p, delta=delta, eta=eta)
        E_frac_avg += E_frac
        E_frac_var += E_frac**2.0
        E_frac_error_avg += E_frac_error
        E_frac_error_var += E_frac_error**2.0

#print(bs*length_scale()/au)        
#Convert to SI
E_ini *= mass_scale()*(length_scale()/time_scale())**2.0
E_thr *= mass_scale()*(length_scale()/time_scale())**2.0
E_imp *= mass_scale()*(length_scale()/time_scale())**2.0
        
#Normalise
E_frac_avg /= N_enc
E_frac_var /= N_enc
E_frac_error_avg /= N_enc
E_frac_error_var /= N_enc

#Calculate variances
E_frac_var -= E_frac_avg**2.0
E_frac_error_var -= E_frac_error_avg**2.0

print('Average fractional energy change =', E_frac_avg)
print('Error on mean of fractional energy change =', E_frac_var**0.5/(N_enc-1)**0.5)
print('Average fractional error on energy change =', E_frac_error_avg)
print('Error on mean of fractional error on energy change =', E_frac_error_var**0.5/(N_enc-1)**0.5)

print('Saving data')
np.savez('impulse_nbody_energy_changes_b10e2_5au_Nenc10e5_record_b_2.npz', E_ini=E_ini, E_thr=E_thr, E_imp=E_imp, b=bs)
print('Finished')
