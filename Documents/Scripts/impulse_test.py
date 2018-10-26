#To test the impulse approximation

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import G, au, parsec, giga, year

from impulse_test_encounter import encounterGrid
from encounters import calc_b_max
from random_direction import randomDirection
from random_binary import setupRandomBinary

from internal_units import *


#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0 / mass_scale()
m2 = 2.0*10.0**30.0 / mass_scale()
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Convert to internal units
rho = rho * (length_scale())**3.0/mass_scale()
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0 / mass_scale()
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0 / length_scale() * time_scale()
#Number density of perturbers
n_p = rho/M_p

#Semi-major axes
a_min = 10.0**3.0 * au / length_scale()
a_max = 10.0**5.0 * au / length_scale()
#Number of a's to test
N_a = 10
#Impact parameters
#b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
#b_max = calc_b_max(M_p, v_rms, a_max, m1, m2)
b_min = a_min
b_max = a_max
#Number of b's to test
N_b = 10

#Number of encounters per each pair of values
#TAKES 2000 MINUTES TO RUN for 20, 20, 1000
N_enc = 10



dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])
a_frac_avg, E_frac_avg, a_bins, b_bins = encounterGrid(m1, m2, v_rms, e, M_p, a_min, a_max, N_a, b_min, b_max, N_b, N_enc)

#Save
np.savez('impulsecontour_hypenc_delta.npz', a_frac_avg=a_frac_avg, E_frac_avg=E_frac_avg, a_bins=a_bins, b_bins=b_bins)
'''
#Load
a_frac_avg = np.load('impulsecontour_2BHTenc_delta.npz')['a_frac_avg']
E_frac_avg = np.load('impulsecontour_2BHTenc_delta.npz')['E_frac_avg']
a_bins = np.load('impulsecontour_2BHTenc_delta.npz')['a_bins']
b_bins = np.load('impulsecontour_2BHTenc_delta.npz')['b_bins']

#Contour plot
#Linear
plt.title('Average fractional error in semi-major axis due to impulse approximation (hyperbolic equations)', wrap=True)
#plt.contourf(a_bins/au, b_bins/au, np.transpose(a_frac_avg))
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), vmin=0.0, vmax=0.06)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
#plt.colorbar()
plt.show()
'''

'''
#Symlog
plt.title('Average fractional error in semi-major axis due to impulse approximation')
ax = plt.gca()
lt = np.min(np.absolute(a_frac_avg[np.nonzero(a_frac_avg)]))
pcm = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(a_frac_avg), norm=colors.SymLogNorm(linthresh=lt))
plt.colorbar(pcm)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

'''
'''
#Log absolute value
plt.title('Absolute average fractional error in semi-major axis due to impulse approximation (double Bahcall et al. equations)', wrap=True)
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.title('Absolute average fractional error in semi-major axis due to impulse approximation (double Bahcall et al. equations)', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), norm=colors.LogNorm())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

#Log absolute value
plt.title('Absolute average fractional error in energy change due to impulse approximation (hyperbolic equations)', wrap=True)
ax = plt.gca()
norm= colors.LogNorm(vmin=10.0**(-7.0),vmax=10.0**2.0)
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(E_frac_avg)), norm=norm)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
plt.title('Absolute average fractional error in energy change due to impulse approximation (hyperbolic equations)', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins*length_scale()/au, b_bins*length_scale()/au, np.transpose(np.absolute(E_frac_avg)), norm=colors.LogNorm(), vmin=0.001, vmax=1000.0)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

'''
plt.title('Sign of average fractional error in semi-major axis due to impulse approximation')
plt.contourf(a_bins/au, b_bins/au, np.transpose(np.sign(a_frac_avg)))
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()

'''
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


