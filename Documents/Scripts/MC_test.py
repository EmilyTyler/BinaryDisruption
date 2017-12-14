#Evolves single binary over time using MC method
#Evolution of the semi-major axis of a binary star due encounters with primordial balck holes (PBHs) which are assumed to have a constant uniform number density

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
import matplotlib.pyplot as plt

#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0**10.0*365.25*24.0*60.0*60.0
#Number of timesteps
N_t = 100
#Time array
t = np.linspace(0.0, T, N_t)
dt = t[1] - t[0]

#Semi-major axis array
A = np.zeros(N_t)
A[0] = 0.1 * 3.086*10.0**16.0
#Eccentricity array
es = np.zeros(N_t)
es[0] = 0.7
#Implement encounters
#The MCEncounters function has been written to move a distribution of semi-major axes forward in time to see what the distribution of semi-major axes now can tell us about the mass of the PBHs. Here I am testing it on one binary and looking at how the semi-major axis changes over time to make sure it makes sense. 
for i in range(1, N_t):
        ([A[i]], [es[i]]) = MCEncounters(v_rms, n_p, dt, m1, m2, M_p, np.array([A[i-1]]), np.array([es[i-1]]), 1)
        
#Plot semi-major axis against time
plt.title('Evolution of semi-major axis over time due to encounters with PBHs')
plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0),A/(3.086*10.0**16.0))
plt.xlabel('Time/Myr')
plt.ylabel('Semimajor axis/pc')
plt.show()
