#To compare my simulations with WSW B19

import numpy as np
import os
os.system("python setup.py build_ext --inplace")
import matplotlib.pyplot as plt

from monte_carlo import MCEncounters_new, MCEncounters_old, WSWEncounters, JTEncounters
from encounters import calc_b_max
from scipy.constants import parsec, au, giga, year, G
from scipy.special import erf

#Mass of binary
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
M_b = m1+m2
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#Relative velocity dispersion
v_rel = np.sqrt(2.0/3.0) * 1.0*10.0**5.0
#Density of dark matter halo solar masses/pc**3
rho = 0.1
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Number density of perturbers
n_p = rho/M_p
#Total simulation time
T = 10.0*giga*year
#Initial semi-major axis
#Semi-major axis separating diffusive and catastrophic regimes
a_crit = G/M_b*(M_p/(0.03*v_rel))**2.0
print('a_crit, au =', a_crit/au)
a_0 = 0.1*parsec
print('a_0, au =', a_0/au)
#Eccentricity
e_0 = 0.7
#Total number of binaries
N_bin = 1000
print('N_bin =', N_bin)
#Number of time steps
N_t = 50
#Time array
dt = T/(N_t-1)
t = np.array([i*dt for i in range(N_t)])
#Number broken
N_broken = np.zeros(N_t)


'''
#Mean time between encounters that will break the binary
#For M_p=1000M_sol breaking impact parameter
b_1 = 1.5*(G*M_p**2.0*a_0**3.0/((m1+m2)*v_rel**2.0))**(1.0/4.0)
#Rate of encounters at b<=b_1
rate_leq = np.pi * b_1**2.0 * n_p * v_rel
#Print mean time between encounters
print('Mean time between encounters that will break the binary, Gyr:', (rate_leq*giga*year)**(-1.0))
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rel, a_0, m1, m2)
#Rate of encounters at b>b_1
rate_g = np.pi * (b_max**2.0-b_1**2.0) * n_p * v_rel
#Energy of binaries
E = G*(m1+m2)/(2.0*a_0)
#Average energy change due to an encounters
delta_E_avg = 7.0/3.0 * (G*M_p*a_0/v_rel)**2.0 * (b_1**(-2.0)-b_max**(-2.0))/(b_max**2.0-b_1**2.0)
print('Mean time for binary to break due to cumulative effects, Gyr:', E/(rate_g*delta_E_avg)/(giga*year))
'''

#Timescales from Binney and Tremaine page 668
print('Perturber mass separating catastrophic and diffusive regimes, M_sol:', 0.03*(v_rel**2.0*M_b*a_0/G)**0.5/(2.0*10.0**30.0))
print('Timescale for catastrophic disruption, Gyr:', 0.07*M_b**0.5/(G**0.5*rho*a_0**1.5)/(giga*year))
print('Timescale for diffusive disruption, Gyr:', 0.002*v_rel*M_b/(G*M_p*rho*a_0)/(giga*year))

a = a_0 + np.zeros(N_bin)
e = e_0 + np.zeros(N_bin)
for i in range(1,N_t):
        print('Timestep',i,'of',N_t-1)
        a, e, N_broken[i] = MCEncounters_new(v_rel, n_p, t[i]-t[i-1], m1, m2, M_p, a, e, np.size(a), a_T=parsec, prefactor=1.0)
        e = e[np.where(a>0.0)]
        a = a[np.where(a>0.0)]
        #a_fin_WSW, e_fin_WSW, N_broken_WSW[i] = WSWEncounters(v_rel, n_p, t[i], m1, m2, M_p, a, e, N_bin)
        #a_fin_JT, e_fin_JT, N_broken_JT[i] = JTEncounters(v_rel, n_p, t[i], m1, m2, M_p, a, e, N_bin)

#Make it cumulative
for i in range(2, N_t):
        N_broken[i] += N_broken[i-1]
#Normalise
N_broken /= N_bin

np.savez('BHTfig2_{}bin_mysim.npz'.format(N_bin), N_broken=N_broken, t=t)

#Find x_0
E_e = G*(m1+m2)/(10.0**3.0*parsec)
x_0 = G*(m1+m2)/(2.0*a_0*E_e)
#Coulomb logarithm ln(b_max/b_FP)
coul_log = 1.0
#e_1/P in Weinberg et al
epsilon = 8.0*(2.0*np.pi)**0.5*n_p*G**2.0*M_p**2.0/v_rel * coul_log
#Dimensionless time
tau = 2.0*epsilon*t/(3.0*E_e)
#Semi-major axis required for x and tau to be similar
print('a for x sim tau, au =', 3.0*G*(m1+m2)/(4.0*epsilon*T)/au)
print('x_0/tau (should be finite) =', x_0/tau[1:])
print('tau (should be very large) =', tau[1:])
#Survival probability
P = np.zeros(N_t)
P[0] = 1.0
P[1:] = erf((x_0/tau[1:])**0.5) - 4.0/(3.0*np.pi**0.5)*(3.0/2.0*(x_0/tau[1:])**0.5 + (x_0/tau[1:])**1.5)*np.exp(-x_0/tau[1:])

#Plot
#plt.plot(t/(giga*year), P, label='Weinberg B19')
#plt.plot(t/(giga*year), 1.0 - N_broken_WSW, label='Weinberg simulation')
#plt.plot(t/(giga*year), 1.0 - N_broken_JT, label='Jiang and Tremaine simulation')
plt.plot(t/(giga*year), 1.0 - N_broken)
#plt.legend()
plt.xlabel('Time, Gyr')
plt.ylabel('Survival probability')
plt.title(r'Predicted survival probability or fraction out of {} binaries that survived as a function of time with $M_p=${}$M_\odot$, $a_\mathrm{{initial}}=${}au'.format(N_bin, int(M_p/(2.0*10.0**30.0)), int(a_0/au)), wrap=True)
plt.show()




























