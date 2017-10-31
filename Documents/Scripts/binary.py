import numpy as np

from matplotlib import pyplot as plt

from encounters import noEncounters
from encounters import binning
from random_binary import setupRandomBinary

#Initialise variables
#Semi-major axis, m
a = 0.1 * 3.086*10.0**16.0
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/Mpc**3
rho = 0.08
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**22.0)**3.0)
#Number of time steps
N_t = 10000
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0

#Number density of perturbers
n_p = rho/M_p
#Global variables
G = 6.67 * 10.0**(-11.0)


#Create array to store positions and velocities over time
X = np.zeros([N_t,4,3])
t = np.zeros(N_t)
#Initialise binary
X[0] = setupRandomBinary(a, e, m1, m2)

#Keep track of a over time
A = np.zeros(N_t)
A[0] = a

#Evolve orbit
#No encounters
#(t, X, A) = noEncounters(N_t, t, X, A, m1, m2)
#Binning method
(t, X, A) = binning(a, v_rms, n_p, N_t, t, X, A, m1, m2, M_p)
#Monte Carlo method



#Plot semi-major axis against time
plt.plot(t,A)
print('A[1]-A[0] = ', (A[1]-A[0]))
print('A[2]-A[1] = ', (A[2]-A[1]))
print('A[3]-A[2] = ', (A[3]-A[2]))
plt.show()
#Plot relative velocity against time
plt.plot(t, np.linalg.norm((X[:,2]-X[:,3]), axis=1))
plt.show()
#Plot relative x position against relative y position
plt.plot((X[:,0,1]-X[:,1,1]), (X[:,0,0]-X[:,1,0]))
plt.show()







