import numpy as np
import math
import random

from matplotlib import pyplot as plt
from scipy.integrate import quad, dblquad

from eccentric_anomaly import findEccentricAnomaly
from evolve_binary import evolveBinary


#Initialise variables
#Semi-major axis, pc
a = 0.1
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/Mpc**3
rho = 0.08
#Number of time steps
N_t = 10000
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, km/s
v_rms = 100.0


#Global variables
G = 6.67 * 10.0**(-11.0)



#Convert to SI
a = a * 3.086*10.0**16.0
v_rms = v_rms * 1000.0
rho = rho * 2.0*10.0**30.0/((3.086*10.0**22.0)**3.0)

#Function to find perturber velocity
def  relativeVelocity():
        #Velocity is in z-direction
        v_rel = np.array([0.0,0.0,v_rms])
        return v_rel


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
#Initial coordinates of second star (cartesian)
x2 = np.array([r*math.cos(f), r*math.sin(f), 0.0])
#Initial velocity of second star
v2 = np.array([- n*a/(math.sqrt(1.0-e**2.0)) * math.sin(f), n*a/(math.sqrt(1.0-e**2.0)) * (e + math.cos(f)), 0.0])
#Initial velocity of first star
v1 = -m2/m1 * v2


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
sini = random.uniform(0.0, 1.0)
i = math.asin(sini)
R_i = np.array([[1.0, 0.0, 0.0],
               [0.0, math.cos(i), math.sin(i)],
               [0.0, -math.sin(i), math.cos(i)]])
x2 = np.transpose(np.dot(R_i, np.transpose(x2)))
v1 = np.transpose(np.dot(R_i, np.transpose(v1)))
v2 = np.transpose(np.dot(R_i, np.transpose(v2)))


#Create array to store positions and velocities over time
X = np.zeros([N_t,4,3])
X[0] = [x1, x2, v1, v2]
t = np.zeros(N_t)



#Temporary time step set up
dt = 0.001 * 2.0*math.pi*math.sqrt(a**3.0/(G*M_b))


#Integrand for collision rates
def integrand(b, v, n_p, v_rms):
        return ((2.0*(2.0*math.pi)**0.5*n_p*b*v**3.0)/(v_rms**3.0)*np.exp(-v**2.0/(2.0*v_rms**2.0)))

#Keep track of a over time
A = np.zeros(N_t)
A[0] = a

for i in xrange(1, N_t):
        
        '''
        #Calculate encounter rates
        #Catastrophic
        k_cat = 0.07
        R_cat = G**0.5*rho*a**(3.0/2.0)/(M_b**(0.5)*k_cat)
        #Diffusive
        Lambda = v_rms**2.0*a/(G*M_p)
        k_diff = 0.022/(math.log(Lambda))
        R_diff = G*M_p*rho*a/(k_diff*v_rms*M_b)
        
        
        #Calculate time step
        dt = 0.1/(np.max([R_cat, R_diff]))
        '''
        #Add time step to time array
        t[i] = t[i-1]+dt
        '''
        
        #Check if there are any encounters
        N_cat = np.random.poisson(R_cat * dt)
        N_diff = np.random.poisson(R_diff * dt)
        print 'N_cat = ', N_cat
        print 'N_diff = ', N_diff
        
        #Implement diffusive encounters
        '''
        
        
        #Evolve orbit
        X[i] = evolveBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
       
        A[i] = G*M_b*np.linalg.norm(X[i,0]-X[i,1])/(2.0*G*M_b-np.linalg.norm(X[i,0]-X[i,1])*np.dot((X[i,2]-X[i,3]),(X[i,2]-X[i,3])))


#Plot position against time
#plt.plot(X[:,0,0],X[:,0,1])
#plt.plot(X[:,1,0],X[:,1,1])
plt.plot(t,A)
print 'A[1]-A[0] = ', (A[1]-A[0])
print 'A[2]-A[1] = ', (A[2]-A[1])
print 'A[3]-A[2] = ', (A[3]-A[2])
plt.show()
plt.plot(t, np.linalg.norm(X[:,2], axis=1))
plt.plot(t, np.linalg.norm(X[:,3], axis=1))
plt.show()
#plt.plot(t, np.linalg.norm(X[:,0]-X[:,1], axis=1))
plt.plot((X[:,0,1]-X[:,1,1]), (X[:,0,0]-X[:,1,0]))
plt.show()







