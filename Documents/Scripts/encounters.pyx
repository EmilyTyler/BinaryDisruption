#Functions to calculate encounter rates and to implement encounters

import numpy as np
from matplotlib import pyplot as plt


from evolve_binary import integrateBinary
from evolve_binary import analyticBinary
from orbital_elements import semimajorAxis
from random_binary import setupRandomBinary
from orbital_elements import orbitalElements
from random_direction import randomDirection

#Global variables
G = 6.67 * 10.0**(-11.0)


#Encounter rate for impact parameters between b0 and b1 and for relative velocities between v0 and v1
def encounterRate(n_p, v_rms, b0, b1, v0, v1):
        rate = np.sqrt(2.0*np.pi)*n_p/v_rms*(b1**2.0-b0**2.0)*((v0**2.0+2.0*v_rms**2.0)*np.exp(-v0**2.0/(2.0*v_rms**2.0))-(v1**2.0+2.0*v_rms**2.0)*np.exp(-v1**2.0/(2.0*v_rms**2.0)))
        return rate


#Evolve binary without encounters
def noEncounters(N_t, t, X, A, m1, m2):
        
        for i in range(1, N_t):
                #Time step
                dt = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Add time step to time array
                t[i] = t[i-1]+dt
                #Evolve orbit
                X[i] = integrateBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)

        return(t, X, A)


#Implements encounters with binning method
def binning(v_rms, n_p, N_t, t, A, es, m1, m2, M_p):
             
        #Mean motion
        n = np.sqrt(G*(m1+m2)/(A[0]**3.0))
        
        #Set up b array
        b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
        #print('b_min = ', b_min)
        b_max = np.max([v_rms/n, 1000.0*b_min])
        #print('b_max = ', b_max)
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
                                
                #Calculate time step
                dt = 1.0/np.amax(R)
                #Add time step to time array
                t[i] = t[i-1] + dt
                
                #Check if there are any encounters
                #Number of encounters matrix:
                N = np.fromfunction(lambda i,j: np.random.poisson(R[i,j]*dt), (N_b,N_v), dtype=int)
                #Array of indices where encounters happen
                i_enc = np.array(np.nonzero(N))
                #Implement encounters
                if np.size(i_enc[0]) > 0:
                        for k in range(np.size(i_enc[0])):
                                for l in range(N[i_enc[0,k], i_enc[1,k]]):
                                        (A[i], es[i]) = encounter(m1, m2, v[i_enc[1,k]], b[i_enc[0,k]], A[i-1], es[i-1], M_p)
                else:
                        A[i] = A[i-1]
                        es[i] = es[i-1]

        return (t, A, es)

#Implement encounters with relative velocity v and impact parameter b using impulse approximation, M_p is perturber mass
def encounter(m1, m2, v, b, a, e, M_p):
        #print('ENCOUNTER!')      
        #Star masses
        m = np.array([m1, m2])
        #90 degree deflection radius
        b_90 = G*(M_p+m)/v**2.0                                                                                                                                                                                                                     
        #Find perturber velocity
        v_vec = v * randomDirection()
        #Open binary
        X = setupRandomBinary(a, e, m1, m2)
        #Centre of mass vector
        R = (m1*X[0] + m2*X[1])/(m1 + m2)
        #Find impact parameter vector
        b_vec = np.dot(R,v_vec)/np.dot(v_vec,v_vec)*v_vec - R
        b_vec = b * b_vec/np.linalg.norm(b_vec)
        #Implement encounter for both stars  
        for i in [0,1]:
                #Calculate impact parameter for this star
                b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/np.dot(v_vec,v_vec) * v_vec + b_vec - X[i]
                #Calculate velocity change in -b direction
                v_perp = 2.0*M_p*v/(m[i]+M_p) * (np.linalg.norm(b_star)/b_90[i])/(1.0 + np.linalg.norm(b_star)**2.0/b_90[i]**2.0) * (-b_star/np.linalg.norm(b_star))
                #Calculate velocity change in -v direction
                v_parr = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + np.linalg.norm(b_star)**2.0/b_90[i]**2.0) * (-v_vec/np.linalg.norm(v_vec))
                #Change velocity
                X[i+2] += v_perp + v_parr                     
        #Close binary
        return orbitalElements(X, m1, m2)









