#To test the impulse approximation
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3

from orbital_elements import orbitalElements
from orbital_elements import semimajorAxis
from random_direction import randomDirection
from random_binary import setupRandomBinary
from evolve_binary import integrateBinary

from scipy.constants import G

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        
        #print('ENCOUNTER!')
        #print("b = ", b)
        #print('V_0 = ', V_0)
        #print('a = ', a)
        #print('e = ', e)
        #Star masses
        m = np.array([m1, m2])                                                                                                                                                                                                                                             
        notBound_imp = False
        a_imp = a
        e_imp = e
        a_frac = 0.0
        e_diff = 0.0
        #If the encounter is not negligible
        if 10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0 > 0.0:             
                #Find perturber velocity
                v_vec = V_0 * randomDirection()
                #Open binary
                X = setupRandomBinary(a, e, m1, m2)
                #Centre of mass vector
                R = (m1*X[0] + m2*X[1])/(m1 + m2)
                #Find impact parameter vector
                b_vec = np.dot(R,v_vec)/V_0**2.0*v_vec - R
                b_vec = b * b_vec/np.linalg.norm(b_vec)
                
                #Impulse approximation:
                #New star velocities for impulse approximation
                V_imp = np.zeros((2,3), dtype=float)               
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
                        #New velocity
                        V_imp[i] = X[i+2] + v_perp + v_parr
        
                #Three body encounter:          
                #Time array
                t = np.array([0.0])
                #Perturber starting distance parameter
                w = np.sqrt(10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0)/V_0
                #End time
                t_end = 2.0*w
                #Initial positions and velocities
                x = np.array([[X[0], X[1], b_vec - w*v_vec, X[2], X[3], v_vec]])
                #Masses
                M = np.array([m1, m2, M_p])
                #Relative separations
                r_12 = np.linalg.norm(x[0,0]-x[0,1])
                r_13 = np.linalg.norm(x[0,0]-x[0,2])
                r_23 = np.linalg.norm(x[0,1]-x[0,2])
                v_1 = np.linalg.norm(x[0,3])
                v_2 = np.linalg.norm(x[0,4])
                v_3 = V_0
                E = np.array([0.5*m1*v_1**2.0+0.5*m2*v_2**2.0+0.5*M_p*v_3**2.0 - G*m1*m2/r_12 - G*m1*M_p/r_13 - G*m2*M_p/r_23])
                delta = np.array([M_p*m[:2]*r_12**2.0/(m1*m2*np.array([r_13**2.0, r_23**2.0]))])
                #Initialise counter
                i = 1
                while t[i-1] < t_end:
                        (x_new, dt) = integrateBinary(3, x[i-1], M)
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
                #Semimajor axis difference
                a_diff = a_imp - a_thr
                a_frac = a_diff/a_imp
                #Eccentricity difference
                e_diff = e_imp - e_thr
                #print('a_frac = ', a_frac)
                #print('e_diff = ', e_diff)
                #Close binary
        return (notBound_imp, a_imp, e_imp, a_frac, e_diff)
        
        
        
        
        
        
        
        
        
        
        
        
        