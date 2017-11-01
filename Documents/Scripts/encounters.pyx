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
def binning(a, v_rms, n_p, N_t, t, X, A, m1, m2, M_p):
             
        #Mean motion
        n = np.sqrt(G*(m1+m2)/(a**3.0))
        
        #Set up b array
        b_min = (np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
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
                
                #Maximum time step
                dt_max = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Calculate time step
                dt = np.min([1.0/np.amax(R),dt_max])
                #Add time step to time array
                t[i] = t[i-1] + dt
                
                #Check if there are any encounters
                #Number of encounters matrix:
                N = np.fromfunction(lambda i,j: np.random.poisson(R[i,j]*dt), (N_b,N_v), dtype=int)
                #
                #if i==50000:
                #        N[9,0] = 1
                #Array of indices where encounters happen
                i_enc = np.array(np.nonzero(N))
                #if i==50000:
                #        print('i_enc = ', i_enc)
                #Implement encounters
                for k in range(np.size(i_enc[0])):
                        for l in range(N[i_enc[0,k], i_enc[1,k]]):
                                X[i-1] = encounter(m1, m2, v[i_enc[1,k]], b[i_enc[0,k]], X[i-1], M_p)                                                
                
                #Evolve orbit
                X[i] = integrateBinary(X[i-1,0], X[i-1,1], X[i-1,2], X[i-1,3], m1, m2, dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)
        return (t, X, A)

#Implement encounters with relative velocity v and impact parameter b using impulse approximation, M_p is perturber mass
def encounter(m1, m2, v, b, X, M_p):
        #print('ENCOUNTER!')      
        #print('b = ', b)
        #print('v = ', v)
        #print('M_p = ', M_p)
        #print('Before: X = ', X)
        #Star masses
        m = np.array([m1, m2])
        #print('m = ', m)
        #90 degree deflection radius
        b_90 = G*(M_p+m)/v**2.0                                                                                                                                                                                                                     
        #print('b_90 = ', b_90)
        #Find which star is closer
        closest_star = np.random.randint(2)
        #print('Closest star = ', closest_star)
        #Calculate impact parameters
        bs = np.zeros(2)
        bs[closest_star] = abs(m[closest_star-1]*(np.linalg.norm(X[0]-X[1]))/(m1+m2) - b)
        bs[closest_star-1] = m[closest_star]*(np.linalg.norm(X[0]-X[1]))/(m1+m2) + b
        #print('bs = ', bs)
        #Implement encounter for both stars  
        for i in [0,1]:
                if 10.0**(-5.0) < bs[i]/b_90[i] < 10.0**5.0:
                        print('b is mid-range')
                v_perp = 2.0*M_p*v/(m[i]+M_p) * (bs[i]/b_90[i])/(1.0 + bs[i]**2.0/b_90[i]**2.0)
                #print('v_perp = ', v_perp)
                v_parr = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + bs[i]**2.0/b_90[i]**2.0)
                #print('v_parr = ', v_parr)
                #Change perpendicular velocity
                if (i == closest_star) and (m[i-1]*np.linalg.norm(X[0]-X[1])/(m1+m2) > b):
                        #print('Perturber between stars')
                        X[i+2] += (X[i]-X[i+1])/np.linalg.norm(X[i]-X[i+1])*v_perp
                else: 
                        X[i+2] += (X[closest_star-1]-X[closest_star])/np.linalg.norm(X[closest_star-1]-X[closest_star])*v_perp
                #Change parallel velocity
                #ASSUME VELOCITY IS IN X-Y PLANE FOR NOW
                X[i+2] += np.cross([0.0, 0.0, 1.0], (X[closest_star-1]-X[closest_star]))/(np.linalg.norm(np.cross([0.0, 0.0, 1.0], (X[closest_star-1]-X[closest_star])))) * (-v_parr)               
        #print('After: X = ', X)      
        return X













