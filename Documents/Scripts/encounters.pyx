#Functions to calculate encounter rates and to implement encounters

import numpy as np
cimport numpy as np 
from matplotlib import pyplot as plt


from evolve_binary import integrateBinary
from evolve_binary import analyticBinary
from orbital_elements import semimajorAxis
from random_binary import setupRandomBinary
from orbital_elements import orbitalElements
from orbital_elements import notBound
from random_direction import randomDirection

# cython profile=True

#Global variables
cdef double G
G = 6.67 * 10.0**(-11.0)


#Encounter rate for impact parameters between b0 and b1 and for relative velocities between v0 and v1
def encounterRate(double n_p, double v_rms, double b0, double b1, double v0, double v1):
        cdef double rate
        rate = np.sqrt(2.0*np.pi)*n_p/v_rms*(b1**2.0-b0**2.0)*((v0**2.0+2.0*v_rms**2.0)*np.exp(-v0**2.0/(2.0*v_rms**2.0))-(v1**2.0+2.0*v_rms**2.0)*np.exp(-v1**2.0/(2.0*v_rms**2.0)))
        return rate


#Evolve binary without encounters
def noEncounters(int N_t, np.ndarray t, np.ndarray X, np.ndarray A, double m1, double m2):
        cdef int i
        cdef double dt
        for i in range(1, N_t):
                #Time step
                dt = 0.0005 * 2.0*np.pi*np.sqrt(A[i-1]**3.0/(G*(m1+m2)))
                #Add time step to time array
                t[i] = t[i-1]+dt
                #Evolve orbit
                X[i] = integrateBinary(2, X[i-1], np.array([m1,m2]), dt)
                #Semi-major axis
                A[i] = semimajorAxis(X[i], m1, m2)
        return(t, X, A)


#Implements encounters with binning method
def binning(double v_rms, double n_p, int N_t, np.ndarray[double, ndim=1] t, np.ndarray[double, ndim=1] A, np.ndarray[double, ndim=1] es, double m1, double m2, double M_p):
        
        #Mean motion
        cdef double n = np.sqrt(G*(m1+m2)/(A[0]**3.0))
        
        #Set up b array
        cdef double b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
        #print('b_min = ', b_min)
        cdef double b_max = np.max([v_rms/n, 1000.0*b_min])
        #print('b_max = ', b_max)
        #Number of bins
        cdef int N_b = 10
        #Width of logarithmically spaced bins
        cdef double dlogb = (np.log(b_max)-np.log(b_min))/N_b
        cdef np.ndarray b = np.zeros([N_b], dtype=float)
        b = np.fromfunction(lambda i: b_min*np.exp(i*dlogb), (N_b,))

        #Set up v array
        cdef double v_min = 0.001 * v_rms
        cdef double v_max = 1000.0 * v_rms
        #Number of bins
        cdef int N_v = 10
        #Width of logarithmically spaced bins
        cdef double dlogv = (np.log(v_max)-np.log(v_min))/N_v
        cdef np.ndarray v = np.zeros([N_v], dtype=float)
        v = np.fromfunction(lambda i: v_min*np.exp(i*dlogv), (N_v,))

        #Matrix of encounter rates
        #R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
        cdef np.ndarray R = np.zeros([N_b,N_v], dtype=float)
        cdef int i,j
        #R = np.fromfunction(lambda i,j: encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv)), (N_b,N_v), dtype=int)
        for i in range(N_b):
                for j in range(N_v):
                        R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv))

        
        cdef double dt
        cdef np.ndarray N = np.zeros([N_b,N_v], dtype=int)
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
                                        (notBound, (A[i], es[i])) = encounter(m1, m2, v[i_enc[1,k]], b[i_enc[0,k]], A[i-1], es[i-1], M_p)
                else:
                        A[i] = A[i-1]
                        es[i] = es[i-1]
                        notBound = False
                if notBound:
                        print('Binary broken!')
                        t[i:] = [t[i] + dt]*len(A[i:])
                        A[i:] = [np.inf]*len(A[i:])
                        es[i:] = [es[i]]*len(A[i:])
                        break
                
        return (t, A, es)

cdef np.ndarray m = np.zeros(2, dtype=float)
cdef np.ndarray b_90 = np.zeros(2, dtype=float)
cdef np.ndarray v_vec = np.zeros(3, dtype=float)
cdef np.ndarray X = np.zeros((4,3), dtype=float)
cdef np.ndarray R = np.zeros(3, dtype=float)
cdef np.ndarray b_vec = np.zeros(3, dtype=float)
cdef int i
cdef np.ndarray b_star = np.zeros(3, dtype=float)
cdef np.ndarray v_perp = np.zeros(3, dtype=float)
cdef np.ndarray v_parr = np.zeros(3, dtype=float)
cdef double v_vec_norm, b_star_norm

#Implement encounters with relative velocity v and impact parameter b using impulse approximation, M_p is perturber mass
def encounter(double m1, double m2, double v, double b, double a, double e, double M_p):
        #print('ENCOUNTER!')      
        #Star masses
        m = np.array([m1, m2])
        #90 degree deflection radius
        b_90 = G*(M_p+m)/v**2.0                                                                                                                                                                                                                     
        #Find perturber velocity
        v_vec = v * randomDirection()
        v_vec_norm = np.linalg.norm(v_vec)
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
                b_star_norm = np.linalg.norm(b_star)
                #Calculate velocity change in -b direction
                v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-b_star/b_star_norm)
                #Calculate velocity change in -v direction
                v_parr = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/v_vec_norm)
                #Change velocity
                X[i+2] += v_perp + v_parr                     
        #Close binary
        return (notBound(X, m1, m2), orbitalElements(X, m1, m2))









