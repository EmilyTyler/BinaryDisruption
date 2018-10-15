#To test the impulse approximation
import math
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3


from orbital_elements import orbitalElements
from orbital_elements import semimajorAxis
from random_direction import randomDirection
from random_binary import setupRandomBinary
from evolve_binary import integrateBinary
from encounters import calc_b_max, impactAndVelocityVectors
from internal_units import *

G = G()
#from scipy.constants import G, parsec, au, giga, year

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        print('ENCOUNTER!')
        print("b = ", b)
        print('V_0 = ', V_0)
        print('a = ', a)
        print('e = ', e)
        
        #Star masses
        cdef np.ndarray m = np.array([m1, m2]) 
        
        #Initialising variables
        notBound_thr = False
        cdef double a_thr = a
        cdef double e_thr = e
        cdef double a_frac = 0.0
        cdef double e_diff = 0.0
        cdef double E_frac = 0.0
        cdef np.ndarray v_vec, X, b_vec, V_imp, b_90, b_star, v_perp, v_parr, t, x, M, x_new, v_BHT
        cdef double b_vec_norm, b_star_norm, w, t_end, dt
        cdef int i
        
        #If the encounter is not negligible   
        if (10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0 > 0.0):
                #Open binary
                X = setupRandomBinary(a, e, m1, m2)
                print('X = ', X)
                #Find impact parameter vector and velocity vector
                b_vec, v_vec = impactAndVelocityVectors(b, V_0)
                print('b_vec =', b_vec)
                print('v_vec =', v_vec)
                #Impulse approximation:
                #New star velocities for impulse approximation
                V_imp = np.zeros((2,3), dtype=float)               
                #90 degree deflection radius
                b_90 = G*(M_p+m)/V_0**2.0 
                print('b_90 =', b_90)
                for i in range(2):
                        print('i =', i)
                        #Calculate impact parameter for this star
                        b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - X[i]
                        print('b_star = ', b_star)
                        b_star_norm = np.sqrt(b_star[0]**2.0+b_star[1]**2.0+b_star[2]**2.0)
                        print('b_star_norm = ', b_star_norm)
                        #Calculate velocity change in b direction
                        v_perp = 2.0*M_p*V_0/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (b_star/b_star_norm)
                        print('v_perp =', v_perp)
                        #Calculate velocity change in -v_vec direction
                        v_parr = 2.0*M_p*V_0/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/V_0)
                        print('v_parr =', v_parr)
                        if b > a:
                                print('b>a')
                                v_BHT = 2.0*G*M_p*a/(b_star_norm**2.0*V_0) * (b_star/b_star_norm)
                        else:
                                print('b<a')
                                v_BHT = 2.0*G*M_p/(b_star_norm*V_0) * (b_star/b_star_norm)
                        print('v_BHT =', v_BHT)
                        #New velocity
                        #V_imp[i] = X[i+2] + v_perp + v_parr
                        V_imp[i] = X[i+2] + v_BHT
                print('V_imp =', V_imp)
                        
                #Three body encounter:          
                #Time array
                t = np.array([0.0])
                #Perturber starting distance parameter
                w = np.sqrt(10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0)/V_0
                #print('w = ', w)
                #End time
                t_end = 2.0*w
                #Initial positions and velocities
                x = np.array([[X[0], X[1], b_vec - w*v_vec, X[2], X[3], v_vec]])
                #print('x[0] = ', x)
                #Masses
                M = np.array([m1, m2, M_p])
                '''
                #Relative separations
                r_12 = np.linalg.norm(x[0,0]-x[0,1])
                r_13 = np.linalg.norm(x[0,0]-x[0,2])
                r_23 = np.linalg.norm(x[0,1]-x[0,2])
                v_1 = np.linalg.norm(x[0,3])
                v_2 = np.linalg.norm(x[0,4])
                v_3 = V_0
                E = np.array([0.5*m1*v_1**2.0+0.5*m2*v_2**2.0+0.5*M_p*v_3**2.0 - G*m1*m2/r_12 - G*m1*M_p/r_13 - G*m2*M_p/r_23])
                delta = np.array([M_p*m[:2]*r_12**2.0/(m1*m2*np.array([r_13**2.0, r_23**2.0]))])
                '''
                #Initialise counter
                i = 1
                #print('t_end =', t_end)
                while t[i-1] < t_end:
                        #print(x[i-1])
                        (x_new, dt) = integrateBinary(3, x[i-1], M, n=1)
                        x = np.append(x, [x_new], axis=0)
                        t = np.append(t, t[i-1]+dt)
                        #print('dt =', dt)
                        '''
                        #Some checks
                        #Relative separations
                        r_12 = np.linalg.norm(x[i,0]-x[i,1])
                        r_13 = np.linalg.norm(x[i,0]-x[i,2])
                        r_23 = np.linalg.norm(x[i,1]-x[i,2])
                        v_1 = np.linalg.norm(x[i,3])
                        v_2 = np.linalg.norm(x[i,4])
                        v_3 = np.linalg.norm(x[i,5])
                        #Total energy
                        E = np.append(E, 0.5*m1*v_1**2.0+0.5*m2*v_2**2.0+0.5*M_p*v_3**2.0 - G*m1*m2/r_12 - G*m1*M_p/r_13 - G*m2*M_p/r_23)
                        #Fractional PBH force strength
                        delta = np.append(delta, [M_p*m[:2]*r_12**2.0/(m1*m2*np.array([r_13**2.0, r_23**2.0]))], axis=0)
                        '''
                        #Increment counter
                        i += 1
                #print(x[i-1])
                #print(t[i-1])
                '''
                #Plot energy against time
                plt.plot(t, E)
                plt.show()
                #Plot fractional PBH force strength against time
                plt.plot(t, delta[:,0])
                plt.plot(t, delta[:,1])
                plt.show()
                #Plot perturber speed against time
                plt.plot(t, np.linalg.norm(x[:,5], axis=1))
                plt.show()  
                #Plot paths of stars
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot(x[:,0,0], x[:,0,1], x[:,0,2])
                ax.plot(x[:,1,0], x[:,1,1], x[:,1,2])
                ax.plot(x[:,2,0], x[:,2,1], x[:,2,2])
                plt.show()
                '''
                '''
                #Generate animation
                base_interval = 10
                fig = plt.figure()
                ax = p3.Axes3D(fig)
                ax.set_xlim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                ax.set_ylim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                ax.set_zlim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                graph = ax.scatter(x[0,:,0], x[0,:,1], x[0,:,2])
                def update(i):
                        ax = p3.Axes3D(fig)
                        ax.set_xlim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                        ax.set_ylim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                        ax.set_zlim3d([-5.0*10.0**16.0, 5.0*10.0**16.0])
                        graph = ax.scatter(x[i,:,0], x[i,:,1], x[i,:,2])
                        anim.event_source.interval = base_interval * (t[i+1]-t[i])/t[1]
                        return(graph)        
                anim = animation.FuncAnimation(fig, update, interval=base_interval, repeat=True)
                plt.show()
                '''
                
                #Set new velocity
                X[2:] = V_imp 
                print('X =', X)
                
                #New Semi-major axis and eccentricities
                (notBound_imp, a_imp, e_imp) = orbitalElements(X, m1, m2)
                (notBound_thr, a_thr, e_thr) = orbitalElements(np.array([x[i-1,0], x[i-1,1], x[i-1,3], x[i-1,4]]), m1, m2)
                print('a_imp = ', a_imp)
                print('a_thr = ', a_thr)
                
                #Semimajor axis difference
                a_diff = a_imp - a_thr
                a_frac = a_diff/a_thr
                #Eccentricity difference
                e_diff = e_imp - e_thr
                #print('a_frac = ', a_frac)
                #print('e_diff = ', e_diff)
                
                #Reduced energies
                E_imp = -G*(m1+m2)/(2.0*a_imp)
                E_thr = -G*(m1+m2)/(2.0*a_thr)
                print('E_imp =', E_imp)
                print('E_thr =', E_thr)
                
                #E_frac = (E_imp - E_thr)/E_thr               
                
                E_ini = -G*(m1+m2)/(2.0*a)
                print('E_ini =', E_ini)
                delta_E_imp = E_imp - E_ini
                delta_E_thr = E_thr - E_ini
                E_frac = (delta_E_imp - delta_E_thr)/delta_E_thr
                print('delta_E_imp =', delta_E_imp)
                print('delta_E_thr =', delta_E_thr)
                print('E_frac =', E_frac)
                
        return (notBound_thr, a_thr, e_thr, a_frac, e_diff, E_frac)
        
        
def encounterGrid(double m1, double m2, double v_rms, double e, double M_p, double a_min, double a_max, int N_a, double b_min, double b_max, int N_b, int N_enc):
        #Set up logarithmic a bins
        cdef double dloga = (np.log(a_max)-np.log(a_min))/N_a
        cdef np.ndarray a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
        #Set up logarithmic b bins
        cdef double dlogb = (np.log(b_max)-np.log(b_min))/N_b
        cdef np.ndarray b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])
        #Average fractional difference in a
        cdef np.ndarray a_frac_avg = np.zeros((N_a, N_b), dtype=float)
        #Average fractional difference in energy
        cdef np.ndarray E_frac_avg = np.zeros((N_a, N_b), dtype=float)
        
        cdef double a_imp, e_imp, a_frac, e_diff, E_frac
        for i in range(N_a):
                #print('a: ', i+1, ' out of ', N_a)
                for j in range(N_b):
                        #print('b: ',j+1,' out of ',N_b)
                        for k in range(N_enc):
                                #print('Enc: ', k+1, ' out of ', N_enc)
                                (notBound_imp, a_imp, e_imp, a_frac, e_diff, E_frac) = impulseTestEncounter(m1, m2, v_rms, b_bins[j], a_bins[i], e, M_p)
                                a_frac_avg[i,j] += a_frac
                                E_frac_avg[i,j] += E_frac
        #Normalise a_frac_avg
        a_frac_avg /= N_enc
        #Normalise E_frac_avg
        E_frac_avg /= N_enc
        return a_frac_avg, E_frac_avg, a_bins, b_bins

def impulseTestEquations(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        #Star masses
        cdef np.ndarray m = np.array([m1, m2]) 
        
        #Initialising variables
        notBound_thr = False
        cdef double a_thr = a
        cdef double e_thr = e
        cdef double a_frac = 0.0
        cdef double e_diff = 0.0
        cdef double E_frac = 0.0
        cdef np.ndarray v_vec, X, b_vec, b_90, b_star, v_perp, v_parr, v_BHT, v_hyp
        cdef double b_vec_norm, b_star_norm
        
        v_hyp = np.zeros((2,3))
        v_BHT = np.zeros((2,3))

        
        #Open binary
        X = setupRandomBinary(a, e, m1, m2)
        #print('X = ', X)
        #Find impact parameter vector and velocity vector
        b_vec, v_vec = impactAndVelocityVectors(b, V_0)
                 
        #90 degree deflection radius
        b_90 = G*(M_p+m)/V_0**2.0 
        for i in range(2):
                #Calculate impact parameter for this star
                b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - X[i]
                #print('b_star = ', b_star)
                b_star_norm = np.sqrt(b_star[0]**2.0+b_star[1]**2.0+b_star[2]**2.0)
                #print('b_star_norm = ', b_star_norm)
                #Calculate velocity change in b direction
                v_perp = 2.0*M_p*V_0/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (b_star/b_star_norm)
                #Calculate velocity change in -v_vec direction
                v_parr = 2.0*M_p*V_0/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/V_0)
                v_hyp[i] = v_perp + v_parr
                if b > a:
                        v_BHT[i] = 2.0*G*M_p*a/(b_star_norm**2.0*V_0) * (b_star/b_star_norm)
                else:
                        v_BHT[i] = 2.0*G*M_p/(b_star_norm*V_0) * (b_star/b_star_norm)
                       
        #New Semi-major axis and eccentricities
        (notBound_hyp, a_hyp, e_hyp) = orbitalElements(np.array([X[0], X[1], X[2]+v_hyp[0], X[3]+v_hyp[1]]), m1, m2)
        (notBound_BHT, a_BHT, e_BHT) = orbitalElements(np.array([X[0], X[1], X[2]+v_BHT[0], X[3]+v_BHT[1]]), m1, m2)

                
        #Semimajor axis difference
        a_diff = a_hyp - a_BHT
        a_frac = a_diff/a_BHT
        #Eccentricity difference
        e_diff = e_hyp - e_BHT
        #print('a_frac = ', a_frac)
        #print('e_diff = ', e_diff)
        
        #Reduced energies
        E_hyp = -G*(m1+m2)/(2.0*a_hyp)
        E_BHT = -G*(m1+m2)/(2.0*a_BHT)
                 
        
        E_ini = -G*(m1+m2)/(2.0*a)
        delta_E_hyp = E_hyp - E_ini
        delta_E_BHT = E_BHT - E_ini
        E_frac = (delta_E_hyp - delta_E_BHT)/delta_E_BHT
        
        print('delta_E_hyp =', delta_E_hyp)
        print('delta_E_BHT =', delta_E_BHT)
        print('E_frac =', E_frac)
                
        return (notBound_thr, a_thr, e_thr, a_frac, e_diff, E_frac)

                        
        
        
        
        
        
        
        
        
