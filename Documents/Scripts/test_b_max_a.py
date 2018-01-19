import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from encounters import calc_b_max, integrateEncounter
from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors

from scipy.constants import au, parsec

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Mass of perturber
M_p = 10.0 * 2.0*10.0**30.0
#Perturber number density
n_p = rho/M_p

#Semi-major axes
a_min = 10.0**3.0 * au
a_max = 10.0**12.0 * au
#Number of a's to test
N_a = 50
#Number of encounters per each pair of values
N_enc = 50
#Delta
delta=10.0**(-6.0)

#Set up logarithmic a bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
#Average fractional difference in a
a_frac_avg = np.zeros(N_a, dtype=float)

#b_max array
b_max = np.array([calc_b_max(M_p, v_rms, a_bins[i], m1, m2, delta) for i in range(N_a)])

for i in range(N_a):
        for j in range(N_enc):
                (notBound, a_new, e_new) = integrateEncounter(m1, m2, v_rms, b_max[i], a_bins[i], e, M_p)
                a_frac_avg[i] += (a_new-a_bins[i])/a_bins[i]
#Normalise a_frac_avg
a_frac_avg /= N_enc

#Plot absolute a_frac_avg against a and b_max/a
fig = plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
plt.title(r'Absolute average fractional change in $a$ due to a single encounter at $b=b_{{max}}$ with $M_p = {}$, $N_a = {}$, $N_{{enc}} = {}$, log10($\delta$) = {}'.format(M_p/(2.0*10.0**30.0), N_a, N_enc, np.log10(delta)))
ax1.loglog(a_bins/au, np.absolute(a_frac_avg), color='dodgerblue')
ax1.set_xlabel(r'$a$, au')
ax1.set_ylabel(r'Absolute average fractional change in $a$')
ax2.loglog(a_bins/au, b_max/a_bins, color='darkorange')
ax2.set_ylabel(r'$b_\mathrm{max}/a$')
plt.show()

#Plot a_frac_avg against a and b_max/a
fig = plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
plt.title(r'Average fractional change in $a$ due to a single encounter at $b=b_{{max}}$ with $M_p = {}$, $N_a = {{enc}}$, $N_enc = {}$, log10($\delta$) = {}'.format(M_p/(2.0*10.0**30.0), N_a, N_enc, np.log10(delta)))
ax1.plot(a_bins/au, a_frac_avg, color='dodgerblue')
ax1.set_xlabel(r'$a$, au')
ax1.set_ylabel(r'Average fractional change in $a$')
ax1.set_yscale('symlog')
ax1.set_xscale('log')
ax2.loglog(a_bins/au, b_max/a_bins, color='darkorange')
ax2.set_ylabel(r'$b_\mathrm{max}/a$')
plt.show()

















