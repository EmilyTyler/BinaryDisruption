import numpy as np
import os
os.system("python setup.py build_ext --inplace")
import matplotlib.pyplot as plt

from encounters import calc_b_max, encounterRate
from monte_carlo import MCEncounters_new
from scipy.constants import parsec, au, giga, year, G

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
a_0 = 0.1*parsec
#Eccentricity
e_0 = 0.7
#Tidal radius
a_T = parsec
#Number of binaries per simulation
N_bin = 100
print('N_bin =', N_bin)
#Number of simulations
N_sim = 1
#Starting index in file names
i_start = 4

#Time steps
#Minimum impact parameter
b_min = 0.0
#Maximum maximum impact parameter
b_max_max = calc_b_max(M_p, v_rel, a_T, m1, m2)
#Minimum velocity
v_min = 10.0**(-2.0)*v_rel
#Maximum velocity
v_max = 10.0**2.0*v_rel
#Encounter rate
rate = encounterRate(n_p, v_rel, b_min, b_max_max, v_min, v_max)
#Timestep
dt = 0.5/rate
#Number of timesteps
N_t = int(T/dt) + 1
print('N_t =', N_t)
#Adjust timestep
dt = T/(N_t-1)
#Time array
t = np.array([i*dt for i in range(N_t)])

for k in range(N_sim):
        print('Simulation',k+1,'of',N_sim)
        #Initialise
        a = a_0 + np.zeros(N_bin)
        e = e_0 + np.zeros(N_bin)
        N_broken = np.zeros(N_t)
        #Time array
        t = np.array([i*dt for i in range(N_t)])
        #Run simulations
        for i in range(1,N_t):
                a, e, N_broken[i] = MCEncounters_new(v_rel, n_p, t[i]-t[i-1], m1, m2, M_p, a, e, np.size(a), a_T=a_T, prefactor=1.0)
                e = e[np.where(a>0.0)]
                a = a[np.where(a>0.0)]
        print('Filtering')
        #Filter out zeros
        previous_number_zero = False
        for i in range(N_t-1,0,-1):
                if N_broken[i]<1:
                        if previous_number_zero:
                                N_broken[i] = -1
                        previous_number_zero = True        
                else:
                        previous_number_zero = False
        t = t[np.where(N_broken>-1)]
        N_broken = N_broken[np.where(N_broken>-1)]

        #print(N_broken)
        #print(t/(giga*year))
        if (np.size(np.where(N_broken>1)))>0:
                print('Multiple binaries broken:', N_broken[np.where(N_broken>1)])
        #Make it cumulative
        for i in range(2, np.size(N_broken)):
                N_broken[i] += N_broken[i-1]
        #Normalise
        N_broken /= N_bin
        #Save data
        print('Saving')
        np.savez('BHTfig2_{}bin_mysim_BHTenc_{}.npz'.format(N_bin, i_start+k), N_broken=N_broken, t=t)





















