import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt

from encounters import noEncounters
from encounters import binning
from random_binary import setupRandomBinary

from scipy.constants import G

#Initialise variables
#Semi-major axis, m
a = 0.1 * 3.086*10.0**16.0
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.008
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number of time steps
N_t = 500000000
#Mass of perturbers
M_p = 10.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0

#Number density of perturbers
n_p = rho/M_p


#Create time array
t = np.zeros(N_t)

#Keep track of a over time
A = np.zeros(N_t)
A[0] = a

#Keep track of e over time
es = np.zeros(N_t)
es[0] = e

#Evolve orbit
#No encounters
#X = np.zeros((N_t,4,3), dtype=float)
#X[0] = setupRandomBinary(A[0], es[0], m1, m2)
#(t, X, A) = noEncounters(N_t, t, X, A, m1, m2)
#Binning method
(t, A, es) = binning(v_rms, n_p, N_t, t, A, es, m1, m2, M_p)
#(t, A, es, a_frac, e_diff) = binning(v_rms, n_p, N_t, t, A, es, m1, m2, M_p)
#Monte Carlo method

'''
#Try Yoo, Chaname and Gould's equation 5
b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
delta_v = np.sqrt(16.0*np.pi*G**2.0*rho*M_p*t*np.log(a/b_min)/v_rms)
X = setupRandomBinary(a, e, m1, m2)
R = np.linalg.norm(X[0]-X[1])
V = np.linalg.norm(X[2]-X[3])
delta_a = G*R**2.0*delta_v*(m1+m2)*(2.0*V+delta_v)/((R*V**2.0-2.0*G*(m1+m2))*(R*(V+delta_v)**2.0-2.0*G*(m1+m2)))
delta_a += a
'''

#Plot relative x position against relative y position
#plt.plot((X[:,0,1]-X[:,1,1]), (X[:,0,0]-X[:,1,0]))
#plt.show()

#Plot semi-major axis against time
plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0),A/(3.086*10.0**16.0))
#plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0),delta_a/(3.086*10.0**16.0))
plt.xlabel('Time/Myr')
plt.ylabel('Semimajor axis/pc')
plt.show()

#Plot relative velocity against time
#plt.plot(t, np.linalg.norm((X[:,2]-X[:,3]), axis=1))
#plt.show()

#Plot eccentricity over time
plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0), es)
plt.xlabel('Time/Myr')
plt.ylabel('e')
plt.show()

'''
#Plot fractional semimajor axis difference
plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0), a_frac)
plt.xlabel('Time/Myr')
plt.ylabel('Fractional semi-major axis difference')
plt.show()

#Plot eccentricity difference
plt.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0), e_diff)
plt.xlabel('Time/Myr')
plt.ylabel('Eccentricity difference')
plt.show()

#Plot semimajor axis and fractional semimajor axis difference
fig = plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
ax1.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0),A/(3.086*10.0**16.0), color='green')
ax1.set_xlabel('Time/Myr')
ax1.set_ylabel('Semimajor axis/pc')
ax1.tick_params(axis='y', colors='green')
ax2.plot(t/(10.0**6.0*365.25*24.0*60.0*60.0), a_frac, color='indigo')
ax2.set_ylabel('Fractional semi-major axis difference between impulse and three-body integration')
ax2.tick_params(axis='y', colors='indigo')
plt.show()
'''








