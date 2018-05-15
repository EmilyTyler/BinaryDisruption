#To compare my simulations with WSW B19

import numpy as np
import os
os.system("python setup.py build_ext --inplace")
import matplotlib.pyplot as plt

from monte_carlo import MCEncounters
from scipy.constants import parsec, au, giga, year, G
from scipy.special import erf

#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Mass of perturbers
M_p = 10000.0 * 2.0*10.0**30.0
#Relative velocity dispersion
v_rel = 2.2*10.0**5.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Number density of perturbers
n_p = rho/M_p
#Total simulation time
T = 10.0*giga*year
#Initial semi-major axis
a_0 = 10.0**4.0 * au
#Eccentricity
e_0 = 0.7
#Total number of binaries for each time step
N_bin = 1000
#Number of time steps
N_t = 10
#Time array
dt = T/N_t
t = np.array([i*dt for i in range(N_t)])
#Number broken
N_broken = np.zeros(N_t)
for i in range(N_t):
        a = a_0 + np.zeros(N_bin)
        e = e_0 + np.zeros(N_bin)
        a_fin, e_fin, N_broken[i] = MCEncounters(v_rel, n_p, t[i], m1, m2, M_p, a, e, N_bin)
#Normalise
N_broken /= N_bin

#Find x_0
E_e = G*(m1+m2)/(2.0*parsec)
x_0 = G*(m1+m2)/(2.0*a_0*E_e)
#Coulomb logarithm ln(b_max/b_FP)
coul_log = 1.0
#e_1/P in Weinberg et al
epsilon = 4.0*(2.0*np.pi)**0.5*n_p*G**2.0*M_p**2.0/v_rel * coul_log
#Dimensionless time
tau = 2.0*epsilon*t/(3.0*E_e)
#Survival probability
P = np.zeros(N_t)
P[0] = 1.0
P[1:] = erf((x_0/tau[1:])**0.5) - 4.0/(3.0*np.pi**0.5)*(3.0/2.0*(x_0/tau[1:])**0.5 + (x_0/tau[1:])**1.5)*np.exp(-x_0/tau[1:])

#Plot
plt.plot(t/(giga*year), P, label='Weinberg B19')
plt.plot(t/(giga*year), 1.0 - N_broken, label='Simulation')
plt.legend()
plt.xlabel('Time, Gyr')
plt.ylabel('Survival probability')
plt.title(r'Predicted survival probability or fraction out of {} binaries that survived as a function of time with $M_p=${}$M_\odot$, $a_\mathrm{{initial}}=${}au'.format(N_bin, int(M_p/(2.0*10.0**30.0)), int(a_0/au)), wrap=True)
plt.show()




























