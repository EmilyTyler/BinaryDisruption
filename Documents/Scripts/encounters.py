#Functions to calculate encounter rates and to implement encounters

import numpy as np
from scipy.integrate import quad, dblquad
from evolve_binary import evolveBinary

#Global variables
G = 6.67 * 10.0**(-11.0)

#Integrand for encounter rates
def integrand(b, v, n_p, v_rms):
        return ((2.0*(2.0*math.pi)**0.5*n_p*b*v**3.0)/(v_rms**3.0)*np.exp(-v**2.0/(2.0*v_rms**2.0)))


#Encounter rate for impact parameters between b0 and b1 and for relative velocities between v0 and v1
def encounterRate(n_p, v_rms, b0, b1, v0, v1):
        return dblquad(lambda b,v: integrand(b,v,n_p,v_rms), b0, b1, lambda v: v0, lambda v: v1)


#Evolve binary without encounters
def noEncounters(dt, N_t, t, X, A, m1, m2):
        for i in xrange(1, N_t):
                #Add time step to time array
                t[i] = t[i-1]+dt
                #Evolve orbit
                X[i] = evolveBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = G*(m1+m2)*np.linalg.norm(X[i,0]-X[i,1])/(2.0*G*(m1+m2)-np.linalg.norm(X[i,0]-X[i,1])*np.dot((X[i,2]-X[i,3]),(X[i,2]-X[i,3])))
        return(t, X, A)

#Implements encounters with binning method
def binning(a, v_rms, n, n_p, dt_max, N_t, t):
             
        #Set up b array
        b_min = a
        b_max = v_rms/n
        #Number of bins
        N_b = 10
        #Width of logarithmically spaced bins
        dlogb = (np.log(b_max)-np.log(b_min))/N_b
        b = np.fromfunction(lambda i: b_min*np.exp(i*dlogb), (N_b,))

        #Set up v array
        v_min = 0.001 * v_rms
        v_max = 1000.0 * v_rms
        #Number of bins
        N_v = 10
        #Width of logarithmically spaced bins
        dlogv = (np.log(v_max)-np.log(v_min))/N_v
        v = np.fromfunction(lambda i: v_min*np.exp(i*dlogv), (N_v,))

        #Matrix of encounter rates
        #R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
        R = np.fromfunction(lambda i,j: encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv)), (N_b,N_v))
        
        #Calculate time step
        dt = np.min([0.1/np.amax(R),dt_max])
                
        #Add time steps to time array
        t = np.fromfunction(lambda i: i*dt, (N_t,))

        for i in xrange(1, N_t):
                
                #Check if there are any encounters
                #Number of encounters matrix:
                N = np.fromfunction(lambda i,j: np.random.poisson(R[i,j]*dt), (N_b,N_v))
                #Total number of encounters
                N_enc = np.sum(R)
                
                
                
                #Implement encounters
                
                
                
                #Evolve orbit
                X[i] = evolveBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = G*M_b*np.linalg.norm(X[i,0]-X[i,1])/(2.0*G*M_b-np.linalg.norm(X[i,0]-X[i,1])*np.dot((X[i,2]-X[i,3]),(X[i,2]-X[i,3])))

















