#To test the impulse approximation
import math
import numpy as np

from orbital_elements import orbitalElements
from random_direction import randomDirection
from random_binary import setupRandomBinary
from evolve_binary import integrateBinary


#Global variables
G = 6.67 * 10.0**(-11.0)

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        #print('ENCOUNTER!')
        #Star masses
        m = np.array([m1, m2])                                                                                                                                                                                                                    
        #Find perturber velocity
        v_vec = V_0 * randomDirection()
        #Open binary
        X = setupRandomBinary(a, e, m1, m2)
        #Centre of mass vector
        R = (m1*X[0] + m2*X[1])/(m1 + m2)
        #Find impact parameter vector
        b_vec = np.dot(R,v_vec)/V_0**2.0*v_vec - R
        b_vec = b * b_vec/np.linalg.norm(b_vec)
        
        #New star velocities for impulse approximation
        V_imp = np.zeros((2,3), dtype=float)
        #Impulse approximation:
        #90 degree deflection radius
        b_90 = G*(M_p+m)/V_0**2.0 
        for i in range(2):
                #Calculate impact parameter for this star
                b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - X[i]
                b_star_norm = np.linalg.norm(b_star)
                #Calculate velocity change in -b direction
                v_perp = 2.0*M_p*V_0/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-b_star/b_star_norm)
                #Calculate velocity change in -v_vec direction
                v_parr = 2.0*M_p*V_0/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/V_0)
                #Change velocity
                V_imp[i] = X[i+2] + v_perp + v_parr
                      
        #Velocity difference
        cdef double V_diff = 0.0;
        #Three body encounter:       
        if 10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0 > 0.0:
                #Time step
                dt = 0.00005 * 2.0*np.pi*np.sqrt(a**3.0/(G*(m1+m2)))
                #Perturber starting distance parameter
                w = np.sqrt(10.0**4.0*M_p*a**2.0/(np.min(m)) - b**2.0)/(dt*V_0)
                #Number of timesteps
                N_t = int(math.ceil(20*w))
                #print('N_t = ', N_t)
                #Positions and velocities
                x = np.zeros((N_t, 6, 3), dtype=float)
                #Initial positions and velocities
                x[0] = np.array([X[0], X[1], b_vec - 10.0*w*dt*v_vec, X[2], X[3], v_vec])
                #Masses
                M = np.array([m1, m2, M_p])
                for i in range(1, N_t):
                        x[i] = integrateBinary(3, x[i-1], M, dt)
                
                V_thr = np.linalg.norm(x[N_t-1,3] - x[N_t-1,4])
        
                #Velocity difference
                V_diff = np.linalg.norm(V_imp[0] - V_imp[1]) - V_thr

        #else:
                #print('Impulse only')

        X[2:] = V_imp
        #Close binary
        return (orbitalElements(X, m1, m2), V_diff)
        
        
        
        
        
        
        
        
        
        
        
        
        