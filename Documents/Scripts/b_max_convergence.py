#Convergence tests for b_max

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from scipy.constants import au, parsec, giga, year
from monte_carlo import draw_b, draw_vmaxwellian
from encounters import encounterRate, encounter

#Initialise variables
#Semi-major axis, m
a = 10.0**4.0 * au
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = 10.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rel = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#b_max
b_max = 0.1*parsec * (M_p/(2.0*10.0**30.0))**0.5 * ((m1+m2)/(4.0*10.0**30.0))**(-0.25) * (a/(10.0**4.0*au))**0.75 * (v_rel/(2.2*10.0**5.0))**(-0.5)

#Plot delta a against t for a small population
#Number of binaries
N_bin = 10
#Number of time steps
N_t = 10
#Simulation run time
T = 10.0 * giga * year
#Time array
t = np.linspace(0.0, T, N_t)
dt = t[1] - t[0]
#b_max values
b_maxs = np.array([0.1*b_max, b_max, 10.0*b_max])
b_maxs_size = np.size(b_maxs)
b_max_max = np.max(b_maxs)
#Array of a's: A[i,j,k] is the semi-major axis for binary i at time t[j] with b_max=b_maxs[k]
A = np.zeros((N_bin, N_t, b_maxs_size))
#Initialise A
A[:,0,:] += a
#Array of es
es = np.zeros((N_bin, N_t, b_maxs_size))
#Initialise es
es[:,0,:] += e
for i in range(N_bin):
        for j in range(N_t):
                #Calculate number of encounters
                N_enc = np.random.poisson(dt*encounterRate(n_p, v_rel, 0.0, b_max_max, 0.01*v_rel, 100.0*v_rel))
                #Generate b values
                b = draw_b(b_max_max, N_enc)
                #Generate v values
                v = draw_vmaxwellian(v_rel, 0.01*v_rel, 100.0*v_rel, N_enc)
                #Implement encounters for different values of b_max
                for k in range(b_maxs_size):
                        b_subset = b[np.where(b <= b_maxs[k])]
                        v_subset = v[np.where(b <= b_maxs[k])]
                        for l in range(np.size(b_subset)):
                                notBound, A[i,j,k], es[i,j,k] = encounter(m1, m2, v_subset[l], b_subset[l], A[i,j,k], es[i,j,k], M_p)
                                if notBound:
                                        A[i,j,k] = -1.0
                                        es[i,j,k] = -1.0
                                        break
#Find change in semi-major axis
delta_a = np.zeros((N_bin, N_t, b_maxs_size))
delta_a[0] = A[0] - a
for j in range(1,N_t):
        delta_a[j] = A[j] - A[j-1]
#Plot delta_a against time for different b_maxs
for i in range(N_bin):
        for k in range(b_maxs_size):
                
                                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
        
        

