#Functions to calculate encounter rates and to implement encounters

import numpy as np
from evolve_binary import integrateBinary
from evolve_binary import evolveBinary
from orbital_elements import semimajorAxis

#Global variables
G = 6.67 * 10.0**(-11.0)


#Encounter rate for impact parameters between b0 and b1 and for relative velocities between v0 and v1
def encounterRate(n_p, v_rms, b0, b1, v0, v1):
        rate = np.sqrt(2.0*np.pi)*n_p/v_rms*(b1**2.0-b0**2.0)*((v0**2.0+2.0*v_rms**2.0)*np.exp(-v0**2.0/(2.0*v_rms**2.0))-(v1**2.0+2.0*v_rms**2.0)*np.exp(-v1**2.0/(2.0*v_rms**2.0)))
        return rate


#Evolve binary without encounters
def noEncounters(N_t, t, X, A, m1, m2):
        
        for i in range(1, N_t//2):
                #Time step
                dt = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Add time step to time array
                t[i] = t[i-1]+dt
                #Evolve orbit
                X[i] = integrateBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)
        
        for i in range(N_t//2, N_t):
                #Time step
                dt = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Add time step to time array
                t[i] = t[i-1]+dt
                #Evolve orbit
                X[i] = evolveBinary(X[i-1], m1, m2, dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)               
        return(t, X, A)

#Implements encounters with binning method
def binning(a, v_rms, n, n_p, N_t, t, X, A, m1, m2):
             
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
        R = np.fromfunction(lambda i,j: encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv)), (N_b,N_v), dtype=int)


        for i in range(1, N_t):
                
                #Maximum time step
                dt_max = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Calculate time step
                dt = np.min([1.0/np.amax(R),dt_max])
                #Add time step to time array
                t[i] = t[i-1] + dt
                
                #Check if there are any encounters
                #Number of encounters matrix:
                N = np.fromfunction(lambda i,j: np.random.poisson(R[i,j]*dt), (N_b,N_v), dtype=int)
                #Array of indices where encounters happen
                i_enc = np.nonzero(N)
                #Implement encounters
                for k in range(np.size(i_enc[0])):
                        for l in range(N(i[k])):
                                X = encounter(m1, m2, v[i_enc[1,k]], b[i_enc[0,k]], X)                                                
                
                #Evolve orbit
                X[i] = integrateBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)
        return (t, X, A)

#Implement encounters with relative velocity v and impact parameter b using impulse approximation
def encounter(m1, m2, v, b, X):
        print('ENCOUNTER!')
        #IMPLEMENT ENCOUNTER
        return X













