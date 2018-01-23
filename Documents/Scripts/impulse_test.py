#To test the impulse approximation

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import G

from impulse_test_encounter import encounterGrid
from encounters import calc_b_max

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Mass of perturbers
M_p = 100.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Semi-major axes
a_min = 10.0**3.0 * 1.5*10.0**11.0
a_max = 10.0**12.0 * 1.5*10.0**11.0
#Number of a's to test
N_a = 20
#Impact parameters
b_min = (np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
b_max = calc_b_max(M_p, v_rms, a_max, m1, m2)
#Number of b's to test
N_b = 20

#Number of encounters per each pair of values
#TAKES 5.5 HOURS TO RUN for 20, 20, 100
N_enc = 100

a_frac_avg, a_bins, b_bins = encounterGrid(m1, m2, v_rms, e, M_p, a_min, a_max, N_a, b_min, b_max, N_b, N_enc)

#Contour plot
'''
#Linear
plt.title('Average fractional error in semi-major axis due to impulse approximation')
plt.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(a_frac_avg))
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()
'''


#Symlog
plt.title('Average fractional error in semi-major axis due to impulse approximation')
ax = plt.gca()
lt = np.min(np.absolute(a_frac_avg[np.nonzero(a_frac_avg)]))
pcm = ax.pcolormesh(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(a_frac_avg), norm=colors.SymLogNorm(linthresh=lt))
plt.colorbar(pcm)
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()



#Log absolute value
plt.title('Absolute average fractional error in semi-major axis due to impulse approximation')
ax = plt.gca()
cs = ax.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()


plt.title('Sign of average fractional error in semi-major axis due to impulse approximation')
plt.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(np.sign(a_frac_avg)))
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()



#Plot time ratio test for impulse approximation
t_T = np.zeros((N_a, N_b))
for i in range(N_a):
        t_T[i] = np.array([(1.0 + e)/(2.0*np.pi*v_rms)*np.sqrt(G*(m1 + m2)/a_bins[i])]*N_b)               
plt.title('Time taken for PBH to travel between points of closest approach divided by orbital period')
ax = plt.gca()
cs = ax.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(t_T), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()


#Plot deltaV/V for PBH velocity
V_frac = np.zeros((N_b, N_a))
for i in range(N_b):
        V_frac[i] = np.array([2.0*m1*G/(b_bins[i]*v_rms**2.0)]*N_a)
plt.title('PBH fractional velocity difference after an encounter with one star')
ax = plt.gca()
cs = ax.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), V_frac, locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()





















