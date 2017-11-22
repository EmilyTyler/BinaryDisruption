#To test the impulse approximation
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from orbital_elements import orbitalElements
from random_direction import randomDirection
from random_binary import setupRandomBinary
from evolve_binary import integrateBinary


#Global variables
G = 6.67 * 10.0**(-11.0)

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        print('ENCOUNTER!')
        print('b = ', b)
        print('V_0 = ', V_0)
        print('a = ', a)
        print('e = ', e)
        #Star masses
        m = np.array([m1, m2])                                                                                                                                                                                                                                         
        #Velocity difference
        cdef double V_diff = 0.0;     
        notBound = False
        if 10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0 > 0.0:
                
                #Find perturber velocity
                v_vec = V_0 * randomDirection()
                #Open binary
                X = setupRandomBinary(a, e, m1, m2)
                print('X = ', X)
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
                print('V_imp = ', V_imp)        
                #Three body encounter:          
                #Time step
                dt = 0.005 * 2.0*np.pi*np.sqrt(a**3.0/(G*(m1+m2)))
                #Perturber starting distance parameter
                w = np.sqrt(10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0)/(dt*V_0)
                #Number of timesteps
                N_t = int(math.ceil(2*w))
                #print('N_t = ', N_t)
                #Positions and velocities
                x = np.zeros((N_t, 6, 3), dtype=float)
                #Initial positions and velocities
                x[0] = np.array([X[0], X[1], b_vec - w*dt*v_vec, X[2], X[3], v_vec])
                #Masses
                M = np.array([m1, m2, M_p])
                #Checks
                E = np.zeros(N_t)
                E[0] = -G*m1*m2/(2.0*a) + 0.5*M_p*V_0**2.0
                delta = np.zeros((N_t,2))
                delta[0] = M_p*m[:2]*np.linalg.norm(x[0,0]-x[0,1])**2.0/(m1*m2*np.array([np.linalg.norm(x[0,0]-x[0,2])**2.0, np.linalg.norm(x[0,1]-x[0,2])**2.0]))
                for i in range(1, N_t):
                        x[i] = integrateBinary(3, x[i-1], M, dt)
                        #Some checks
                        #Relative separations
                        r_12 = np.linalg.norm(x[i,0]-x[i,1])
                        r_13 = np.linalg.norm(x[i,0]-x[i,2])
                        r_23 = np.linalg.norm(x[i,1]-x[i,2])
                        v_1 = np.linalg.norm(x[i,3])
                        v_2 = np.linalg.norm(x[i,4])
                        v_3 = np.linalg.norm(x[i,5])
                        #Total energy
                        E[i] = 0.5*m1*v_1**2.0+0.5*m2*v_2**2.0+0.5*M_p*v_3**2.0 - G*m1*m2/r_12**2.0 - G*m1*M_p/r_13**2.0 - G*m2*M_p/r_23**2.0
                        #Fractional PBH force strength
                        delta[i] = M_p*m[:2]*r_12**2.0/(m1*m2*np.array([r_13**2.0, r_23**2.0]))
                
                print('x[N_t-1] = ', x[N_t-1])
                t = [i*dt for i in range(N_t)]
                print('t/P = ', t[N_t-1]/(2.0*np.pi*np.sqrt(a**3.0/(G*(m1+m2)))))
                print('P = ', (2.0*np.pi*np.sqrt(a**3.0/(G*(m1+m2)))))
                plt.plot(t, E)
                plt.show()
                plt.plot(t, delta[:,0])
                plt.plot(t, delta[:,1])
                plt.show()
                plt.plot(t, np.linalg.norm(x[:,5], axis=1))
                plt.show()
                #plt.plot(x[:,0,0], x[:,0,1])
                #plt.plot(x[:,1,0], x[:,1,1])
                #plt.plot(x[:,2,0], x[:,2,1])
                #plt.show()
                
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot(x[:,0,0], x[:,0,1], x[:,0,2])
                ax.plot(x[:,1,0], x[:,1,1], x[:,1,2])
                plt.show()
                
                V_thr = np.linalg.norm(x[N_t-1,3] - x[N_t-1,4])
                print('V_thr = ', V_thr)
                #Velocity difference
                V_diff = np.linalg.norm(V_imp[0] - V_imp[1]) - V_thr
                print('V_diff = ', V_diff)
                X[2:] = V_imp     
                #Close binary
                (notBound, a, e) = orbitalElements(X, m1, m2)
                
        return (notBound, a , e, V_diff)
        
        
        
        
        
        
        
        
        
        
        
        
        