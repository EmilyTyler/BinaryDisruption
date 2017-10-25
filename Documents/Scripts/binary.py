import numpy as np
import math
import random

from matplotlib import pyplot as plt

from eccentric_anomaly import findEccentricAnomaly
from encounters import noEncounters
from encounters import binning


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

#Initialise binary
#Total mass of binary
M_b = m1 + m2
#Randomise mean anomaly
M = random.uniform(0.0, 2.0*math.pi)
#Find eccentric anomaly from Kepler's equation
E = findEccentricAnomaly(e, M)
#Find true anomaly
f = 2.0*math.atan(((1.0+e)/(1.0-e))**0.5 * math.tan(0.5*E))
#Initial separation
r = a*(1.0 - e**2.0)/(1.0 + e*math.cos(f))
#Mean motion
n = math.sqrt(G*M_b/(a**3.0))
#Initial coordinates of first star (cartesian)
x1 = np.array([0.0, 0.0, 0.0])
#Initial velocity of first star
v1 = np.array([0.0, 0.0, 0.0])
#Initial coordinates of second star (cartesian)
x2 = np.array([r*math.cos(f), r*math.sin(f), 0.0])
#Initial velocity of second star
v2 = np.array([- n*a/(math.sqrt(1.0-e**2.0)) * math.sin(f), n*a/(math.sqrt(1.0-e**2.0)) * (e + math.cos(f)), 0.0])

'''

#Randomly orient binary
#Rotation about z axis by angle phi
phi = random.uniform(0.0, 2.0*math.pi)
R_phi = np.array([[math.cos(phi), math.sin(phi), 0.0],
                  [-math.sin(phi), math.cos(phi), 0.0],
                  [0.0, 0.0, 1.0]])
x2 = np.transpose(np.dot(R_phi, np.transpose(x2)))
v1 = np.transpose(np.dot(R_phi, np.transpose(v1)))
v2 = np.transpose(np.dot(R_phi, np.transpose(v2)))
#Rotation about x axis by angle i
sinI = random.uniform(0.0, 1.0)
I = math.asin(sinI)
R_I = np.array([[1.0, 0.0, 0.0],
               [0.0, math.cos(I), math.sin(I)],
               [0.0, -math.sin(I), math.cos(I)]])
x2 = np.transpose(np.dot(R_I, np.transpose(x2)))
v1 = np.transpose(np.dot(R_I, np.transpose(v1)))
v2 = np.transpose(np.dot(R_I, np.transpose(v2)))
'''

#Create array to store positions and velocities over time
X = np.zeros([N_t,4,3])
X[0] = [x1, x2, v1, v2]
t = np.zeros(N_t)


#Keep track of a over time
A = np.zeros(N_t)
A[0] = a

#Evolve orbit
#No encounters
(t, X, A) = noEncounters(N_t, t, X, A, m1, m2)
#Binning method
#(t, X, A) = binning(a, v_rms, n, n_p, N_t, t, X, A, m1, m2)
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







