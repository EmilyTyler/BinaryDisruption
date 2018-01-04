#Functions to calculate encounter rates and to implement encounters

import numpy as np
cimport numpy as np 
from matplotlib import pyplot as plt


from evolve_binary import integrateBinary
from orbital_elements import semimajorAxis
from random_binary import setupRandomBinary
from orbital_elements import orbitalElements
from random_direction import randomDirection
from impulse_test_encounter import impulseTestEncounter

from scipy.constants import G


#Encounter rate for impact parameters between b0 and b1 and for relative velocities between v0 and v1
def encounterRate(double n_p, double v_rms, double b0, double b1, double v0, double v1):
        cdef double rate
        rate = np.sqrt(2.0*np.pi)*n_p/v_rms*(b1**2.0-b0**2.0)*((v0**2.0+2.0*v_rms**2.0)*np.exp(-v0**2.0/(2.0*v_rms**2.0))-(v1**2.0+2.0*v_rms**2.0)*np.exp(-v1**2.0/(2.0*v_rms**2.0)))
        return rate

#To find b_max
def calc_b_max(M_p, v, a, m1, m2):
        return (2.0*G*M_p/(v*10.0**(-3.0))*np.sqrt(a/(G*(m1+m2))))

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
def binning(double v_rms, double n_p, double t_end, double a_0, double e_0, double m1, double m2, double M_p):
        cdef int i,j,k
        #Set up b array
        cdef double b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
        #print('b_min = ', b_min)
        cdef double b_max = calc_b_max(M_p, v_rms, a_0, m1, m2)
        #print('b_max = ', b_max)
        #Number of bins
        cdef int N_b = 100
        #Width of logarithmically spaced bins
        cdef double dlogb = (np.log(b_max)-np.log(b_min))/N_b
        cdef np.ndarray b = np.array([b_min*np.exp(i*dlogb) for i in range(N_b)])

        #Set up v array
        cdef double v_min = 0.001 * v_rms
        cdef double v_max = 1000.0 * v_rms
        #Number of bins
        cdef int N_v = 100
        #Width of logarithmically spaced bins
        cdef double dlogv = (np.log(v_max)-np.log(v_min))/N_v
        cdef np.ndarray v = np.array([v_min*np.exp(i*dlogv) for i in range(N_v)])

        #Matrix of encounter rates
        #R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
        cdef np.ndarray R = np.zeros([N_b,N_v], dtype=float)
        #R = np.fromfunction(lambda i,j: encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv)), (N_b,N_v), dtype=int)
        for i in range(N_b):
                for j in range(N_v):
                        R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv))
        #Calculate time step
        cdef double dt = 1.0/np.amax(R)
        #Setup time array
        cdef int N_t = int(np.ceil(t_end/dt + 1))
        cdef np.ndarray t = np.array([dt*i for i in range(N_t)])
        #Array of semi-major axes
        cdef np.ndarray A = np.zeros(N_t, dtype=float)
        A[0] = a_0
        #Array of eccentricities
        cdef np.ndarray es = np.zeros(N_t, dtype=float)
        es[0] = e_0
        #Number of encounters matrix:
        cdef np.ndarray N = np.rollaxis(np.array([[np.random.poisson(R[i,j]*dt, size=N_t) for j in range(N_v)] for i in range(N_b)]), 2)
        #Array of indices where encounters happen
        cdef np.ndarray i_enc = np.transpose(np.array(np.nonzero(N)))
        
        #
        #a_frac = np.zeros((N_t), dtype=float)
        #e_diff = np.zeros((N_t), dtype=float)
        cdef int i_old = 0
        for (i,j,k) in i_enc:
                if i != i_old:
                        A[i] = A[i_old]
                        es[i] = es[i_old]
                #Implement encounters
                (notBound, A[i], es[i]) = encounter(m1, m2, v[k], b[j], A[i], es[i], M_p)
                #(notBound, A[i], es[i], a_frac[i], e_diff[i]) = impulseTestEncounter(m1, m2, v[k], b[j], A[i], es[i], M_p)
                '''
                #Test b_max
                if i != i_old:
                        if abs(A[i]-A[i_old])/A[i_old] > 10.0**(-3.0):
                                print('b_max too small')
                                print('b/b_max = ', b_old/b_max)
                                print('a_frac = ', abs(A[i]-A[i_old])/A[i_old])
                        elif abs(A[i]-A[i_old])/A[i_old] < 10.0**(-13.0):
                                print('b_max too big')
                                print('b/b_max = ', b_old/b_max)
                                print('a_frac = ', abs(A[i]-A[i_old])/A[i_old])
                '''                
                if notBound:
                        print('Binary broken!')
                        t[i:] = [t[i] + dt]*len(A[i:])
                        A[i:] = [np.inf]*len(A[i:])
                        es[i:] = [es[i]]*len(A[i:])
                        break
                i_old = i
                b_old = b[j]
        #Fill in zeros
        while np.array(np.where(A == 0.0)).size > 0:
                for i in np.transpose(np.array(np.where(A == 0.0))):
                        A[i] = A[i-1]
        while np.array(np.where(es == 0.0)).size > 0:
                for i in np.transpose(np.array(np.where(es == 0.0))):
                        es[i] = es[i-1]
        return (t, A, es)
        #return (t, A, es, a_frac, e_diff)

cdef np.ndarray m = np.zeros(2, dtype=float)
cdef np.ndarray b_90 = np.zeros(2, dtype=float)
cdef np.ndarray v_vec = np.zeros(3, dtype=float)
cdef np.ndarray X = np.zeros((4,3), dtype=float)
cdef np.ndarray X_0 = np.zeros((4,3), dtype=float)
cdef np.ndarray R = np.zeros(3, dtype=float)
cdef np.ndarray b_vec = np.zeros(3, dtype=float)
cdef int i
cdef np.ndarray b_star = np.zeros(3, dtype=float)
cdef np.ndarray v_perp = np.zeros(3, dtype=float)
cdef np.ndarray v_parr = np.zeros(3, dtype=float)
cdef double v_vec_norm, b_star_norm

#Implement encounters with relative velocity v and impact parameter b using impulse approximation, M_p is perturber mass
#From Binney and Tremaine hyperbolic encounters
def encounter(double m1, double m2, double v, double b, double a, double e, double M_p):
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
        b_vec = np.dot(R,v_vec)/v**2.0*v_vec - R
        b_vec_norm = np.sqrt(b_vec[0]**2.0 + b_vec[1]**2.0 + b_vec[2]**2.0)
        b_vec = b * b_vec/b_vec_norm
        #Implement encounter for both stars  
        for i in range(2):
                #Calculate impact parameter for this star
                b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/v**2.0 * v_vec + b_vec - X[i]
                b_star_norm = np.sqrt(b_star[0]**2.0 + b_star[1]**2.0 + b_star[2]**2.0)
                #Calculate velocity change in -b direction
                v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-b_star/b_star_norm)
                #print('v_perp = ', v_perp)
                #Calculate velocity change in -v direction
                v_parr = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/v)
                #print('v_parr = ', v_parr)
                #Change velocity
                X[i+2] += v_perp + v_parr                    
        #Close binary
        return orbitalElements(X, m1, m2)


        
        









