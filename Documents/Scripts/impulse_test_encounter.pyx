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
from encounters import calc_b_max

from scipy.constants import G, parsec, au, giga, year

#Function to give the value of b for borderline impulse validity
def bImpulseValid(double a, double n_p, double v_rms, double M_p, double m1, double m2):
        cdef double m = 1.0
        cdef double c = 2.0 - 3.0*m
        cdef double b = parsec * 10.0**(m*np.log10(a/au) + c)
        cdef double b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
        cdef double b_max = calc_b_max(M_p, v_rms, a, m1, m2)
        cdef double answer = b
        if b < b_min:
                answer = b_min
        elif b > b_max:
                answer = b_max
        return answer

def MImpulseValid(a):
        m = 1.0
        c = 1.0 - 9.0*m
        return 2.0*10.0**30.0 * 10.0**(m*np.log10(a/au) + c)

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        #print('ENCOUNTER!')
        #print("b = ", b)
        #print('V_0 = ', V_0)
        #print('a = ', a)
        #print('e = ', e)
        #Star masses
        cdef np.ndarray m = np.array([m1, m2])                                                                                                                                                                                                                                             
        notBound_imp = False
        cdef double a_imp = a
        cdef double e_imp = e
        cdef double a_frac = 0.0
        cdef double e_diff = 0.0

        cdef np.ndarray v_vec, X, R, b_vec, V_imp, b_90, b_star, v_perp, v_parr, t, x, M, x_new
        cdef double b_vec_norm, b_star_norm, w, t_end, dt
        cdef int i
        #If the encounter is not negligible
        if 10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0 > 0.0:             
                #Find perturber velocity
                v_vec = V_0 * randomDirection()
                print('v_vec = ', v_vec)
                #Open binary
                X = setupRandomBinary(a, e, m1, m2)
                print('X = ', X)
                #Centre of mass vector
                R = (m1*X[0] + m2*X[1])/(m1 + m2)
                #Find impact parameter vector
                b_vec = np.dot(R,v_vec)/V_0**2.0*v_vec - R
                b_vec_norm = np.sqrt(b_vec[0]**2.0+b_vec[1]**2.0+b_vec[2]**2.0)
                b_vec = b * b_vec/b_vec_norm
                print('b_vec = ', b_vec)
                #Impulse approximation:
                #New star velocities for impulse approximation
                V_imp = np.zeros((2,3), dtype=float)               
                #90 degree deflection radius
                b_90 = G*(M_p+m)/V_0**2.0 
                for i in range(2):
                        #Calculate impact parameter for this star
                        b_star = (np.dot(X[i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - X[i]
                        print('b_star = ', b_star)
                        b_star_norm = np.sqrt(b_star[0]**2.0+b_star[1]**2.0+b_star[2]**2.0)
                        print('b_star_norm = ', b_star_norm)
                        #Calculate velocity change in b direction
                        v_perp = 2.0*M_p*V_0/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (b_star/b_star_norm)
                        #Calculate velocity change in -v_vec direction
                        v_parr = 2.0*M_p*V_0/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/V_0)
                        #New velocity
                        V_imp[i] = X[i+2] + v_perp + v_parr
                print('V_imp = ', V_imp)
                #Three body encounter:          
                #Time array
                t = np.array([0.0])
                #Perturber starting distance parameter
                w = np.sqrt(10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0)/V_0
                #End time
                t_end = 2.0*w
                #Initial positions and velocities
                x = np.array([[X[0], X[1], b_vec - w*v_vec, X[2], X[3], v_vec]])
                print('x[0] = ', x)
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
                while t[i-1] < t_end:
                        (x_new, dt) = integrateBinary(3, x[i-1], M, n=1)
                        x = np.append(x, [x_new], axis=0)
                        t = np.append(t, t[i-1]+dt)
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
                #ax.plot(x[:,2,0], x[:,2,1], x[:,2,2])
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
                #Semi-major axis and eccentricity differences
                (notBound_imp, a_imp, e_imp) = orbitalElements(X, m1, m2)
                (notBound_thr, a_thr, e_thr) = orbitalElements(np.array([x[i-1,0], x[i-1,1], x[i-1,3], x[i-1,4]]), m1, m2)
                print('x[i-1] = ', x[i-1])
                print('a_imp = ', a_imp)
                print('a_thr = ', a_thr)
                #Semimajor axis difference
                a_diff = a_imp - a_thr
                a_frac = a_diff/a_thr
                #Eccentricity difference
                e_diff = e_imp - e_thr
                #print('a_frac = ', a_frac)
                #print('e_diff = ', e_diff)
                #Close binary
        return (notBound_thr, a_thr, e_thr, a_frac, e_diff)
        
        
def encounterGrid(double m1, double m2, double v_rms, double e, double M_p, double a_min, double a_max, int N_a, double b_min, double b_max, int N_b, int N_enc):
        #Set up logarithmic a bins
        cdef double dloga = (np.log(a_max)-np.log(a_min))/N_a
        cdef np.ndarray a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
        #Set up logarithmic b bins
        cdef double dlogb = (np.log(b_max)-np.log(b_min))/N_b
        cdef np.ndarray b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])
        #Average fractional difference in a
        cdef np.ndarray a_frac_avg = np.zeros((N_a, N_b), dtype=float)

        cdef double a_imp, e_imp, a_frac, e_diff
        for i in range(N_a):
                #print('a: ', i+1, ' out of ', N_a)
                for j in range(N_b):
                        #print('b: ',j+1,' out of ',N_b)
                        for k in range(N_enc):
                                #print('Enc: ', k+1, ' out of ', N_enc)
                                (notBound_imp, a_imp, e_imp, a_frac, e_diff) = impulseTestEncounter(m1, m2, v_rms, b_bins[j], a_bins[i], e, M_p)
                                a_frac_avg[i,j] += a_frac
        #Normalise a_frac_avg
        a_frac_avg /= N_enc
        return a_frac_avg, a_bins, b_bins
        
        
def encounterGrid_M(double m1, double m2, double v_rms, double rho, double e, double M_p_min, double M_p_max, double a_min, double a_max, int N_M, int N_a, int N_enc):
        #Set up logarithmic a bins
        cdef double dloga = (np.log(a_max)-np.log(a_min))/N_a
        cdef np.ndarray a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
        #Set up logarithmic M_p bins
        cdef double dlogM = (np.log(M_p_max)-np.log(M_p_min))/N_M
        cdef np.ndarray M_p_bins = np.array([M_p_min*np.exp(dlogM*i) for i in range(N_M)])
        #Average fractional difference in a
        cdef np.ndarray a_frac_avg = np.zeros((N_a, N_M), dtype=float)

        cdef double a_imp, e_imp, a_frac, e_diff
        for j in range(N_M):
                n_p = rho/M_p_bins[j]
                for i in range(N_a):
                        for k in range(N_enc):
                                (notBound_imp, a_imp, e_imp, a_frac, e_diff) = impulseTestEncounter(m1, m2, v_rms, bImpulseValid( a_bins[i], n_p, v_rms, M_p_bins[j], m1, m2), a_bins[i], e, M_p_bins[j])
                                a_frac_avg[i,j] += a_frac
        #Normalise a_frac_avg
        a_frac_avg /= N_enc
        return a_frac_avg, a_bins, M_p_bins
        
        
        
        
        
        
        
        
