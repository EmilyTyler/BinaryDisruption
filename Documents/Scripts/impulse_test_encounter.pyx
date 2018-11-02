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
from scipy.constants import parsec, au, giga, year

def impulseTestEncounter(double m1, double m2, double V_0, double b, double a, double e, double M_p):
        '''
        print('ENCOUNTER!')
        print("b = ", b)
        print('V_0 = ', V_0)
        print('a = ', a)
        print('e = ', e)
        '''
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
                #print('X = ', X)
                #Find impact parameter vector and velocity vector
                b_vec, v_vec = impactAndVelocityVectors(b, V_0)
                #print('b_vec =', b_vec)
                #print('v_vec =', v_vec)
                #Perturber starting distance parameter
                w = np.sqrt(10.0**6.0*M_p*a**2.0/(np.min(m)) - b**2.0)/V_0
                #print('w =', w)
                #End time
                t_end = 2.0*w
                #print('t_end =', t_end)
                
                
                #Impulse approximation:            
                #90 degree deflection radius
                b_90 = G*(M_p+m)/V_0**2.0 
                b_star = np.zeros((2,3))
                #print('b_90 =', b_90)
                #Initial evolution of binary
                #Time array
                t_imp = np.array([0.0])
                #Initial positions and velocities
                x_imp = np.array([X])
                #Initialise counter
                j = 1
                while t_imp[j-1] <= w:
                        (x_imp_new, dt) = integrateBinary(2, x_imp[j-1], m, n=100, dt_max=w/100.0)
                        x_imp = np.append(x_imp, [x_imp_new], axis=0)
                        t_imp = np.append(t_imp, t_imp[j-1]+dt)                        
                        #Increment counter
                        j += 1
                for i in range(2):
                        #print('i =', i)
                        #Calculate impact parameter for this star
                        b_star[i] = (np.dot(x_imp[j-1,i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - x_imp[j-1,i]
                        #print('b_star = ', b_star)
                        b_star_norm = np.sqrt(b_star[i,0]**2.0+b_star[i,1]**2.0+b_star[i,2]**2.0)
                        #print('b_star_norm = ', b_star_norm)
                        #Calculate velocity change in b direction
                        v_perp = 2.0*M_p*V_0/(m[i]+M_p) * (b_star_norm/b_90[i])/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (b_star[i]/b_star_norm)
                        #print('v_perp =', v_perp)
                        #Calculate velocity change in -v_vec direction
                        v_parr = 2.0*M_p*V_0/(m[i]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[i]**2.0) * (-v_vec/V_0)
                        #print('v_parr =', v_parr)
                        #New velocity
                        x_imp[j-1,i+2] += v_perp + v_parr
                #print('V_imp =', V_imp)
                while t_imp[j-1] < t_end:
                        (x_imp_new, dt) = integrateBinary(2, x_imp[j-1], m, n=100, dt_max=w/100.0)
                        x_imp = np.append(x_imp, [x_imp_new], axis=0)
                        t_imp = np.append(t_imp, t_imp[j-1]+dt)                        
                        #Increment counter
                        j += 1
                N_t_imp = np.size(t_imp)
                #print('t_imp =', t_imp)
                
                
                #Impulse approximation, Bahcall et al.:            
                b_star = np.zeros((2,3))
                #print('b_90 =', b_90)
                #Initial evolution of binary
                #Time array
                t_BHT = np.array([0.0])
                #Initial positions and velocities
                x_BHT = np.array([X])
                #Initialise counter
                j = 1
                while t_BHT[j-1] <= w:
                        (x_BHT_new, dt) = integrateBinary(2, x_BHT[j-1], m, n=100, dt_max=w/100.0)
                        x_BHT = np.append(x_BHT, [x_BHT_new], axis=0)
                        t_BHT = np.append(t_BHT, t_BHT[j-1]+dt)                        
                        #Increment counter
                        j += 1
                for i in range(2):
                        #print('i =', i)
                        #Calculate impact parameter for this star
                        b_star[i] = (np.dot(x_BHT[j-1,i],v_vec) - np.dot(b_vec,v_vec))/V_0**2.0 * v_vec + b_vec - x_BHT[j-1,i]
                        #print('b_star = ', b_star)
                for i in range(2):
                        b_star_norm = np.sqrt(b_star[i,0]**2.0+b_star[i,1]**2.0+b_star[i,2]**2.0)
                        #print('b_star_norm = ', b_star_norm)
                        if b > a:
                                #print('b>a')
                                v_BHT = 2.0*G*M_p/(b**2.0*V_0) * (np.dot(v_vec,x_BHT[j-1,i])*v_vec/(V_0**2.0)+2.0*np.dot(b_vec,x_BHT[j-1,i])*b_vec/(b**2.0) - x_BHT[j-1,i])
                        else:
                                #print('b<a')
                                v_BHT = 2.0*G*M_p/(b_star_norm*V_0) * (b_star[i]/b_star_norm)
                        print('v_BHT =', v_BHT)
                        print('1st term =', 2.0*G*M_p/(b**2.0*V_0) * np.dot(v_vec,x_BHT[j-1,i])*v_vec/(V_0**2.0))
                        print('2nd term =', 2.0*G*M_p/(b**2.0*V_0) * 2.0*np.dot(b_vec,x_BHT[j-1,i])*b_vec/(b**2.0))
                        print('3st term =', 2.0*G*M_p/(b**2.0*V_0) * - x_BHT[j-1,i])
                        print('BHT equa =', 2.0*G*M_p*a/(b**2.0*V_0) * -x_BHT[j-1,i]/np.linalg.norm(x_BHT[j-1,i]))
                        print('norm     =', np.linalg.norm(v_BHT))
                        print('norm BHT =', G*M_p*a/(b**2.0*V_0))
                        #New velocity
                        x_BHT[j-1,i+2] += v_BHT
                #print('V_imp =', V_imp)
                while t_BHT[j-1] < t_end:
                        (x_BHT_new, dt) = integrateBinary(2, x_BHT[j-1], m, n=100, dt_max=w/100.0)
                        x_BHT = np.append(x_BHT, [x_BHT_new], axis=0)
                        t_BHT = np.append(t_BHT, t_BHT[j-1]+dt)                        
                        #Increment counter
                        j += 1
                N_t_BHT = np.size(t_BHT)
                #print('t_imp =', t_imp)
                        
                        
                #Three body encounter:          
                #Time array
                t_thr = np.array([0.0])
                #Initial positions and velocities
                x = np.array([[X[0], X[1], b_vec - w*v_vec, X[2], X[3], v_vec]])
                #print('x[0] = ', x)
                #Masses
                M = np.array([m1, m2, M_p]) 
                #Initialise counter
                i = 1
                while t_thr[i-1] < t_end:
                        (x_new, dt) = integrateBinary(3, x[i-1], M, n=100)
                        x = np.append(x, [x_new], axis=0)
                        t_thr = np.append(t_thr, t_thr[i-1]+dt)                        
                        #Increment counter
                        i += 1
                N_t_thr = np.size(t_thr)
                #print('t_thr =', t_thr)
                
                '''
                #Move everything into binary centre of mass velocity frame
                for i in range(N_t_thr):
                        v = (m1*x[i,3] + m2*x[i,4])/(m1+m2)
                        x[i,3] -= v
                        x[i,4] -= v
                        x[i,5] -= v
                for i in range(N_t_imp):
                        v = (m1*x_imp[i,2] + m2*x_imp[i,3])/(m1+m2)
                        x_imp[i,2] -= v
                        x_imp[i,3] -= v
                for i in range(N_t_BHT):
                        v = (m1*x_BHT[i,2] + m2*x_BHT[i,3])/(m1+m2)
                        x_BHT[i,2] -= v
                        x_BHT[i,3] -= v
                '''
                
                #Calculations
                #Initial energy of binary
                E_ini = -G*m1*m2/(2.0*a)
                #New Semi-major axis and eccentricities
                (notBound_imp, a_imp, e_imp, E_imp) = orbitalElements(x_imp[N_t_imp-1], m1, m2)
                (notBound_thr, a_thr, e_thr, E_thr) = orbitalElements(np.array([x[N_t_thr-1,0], x[N_t_thr-1,1], x[N_t_thr-1,3], x[N_t_thr-1,4]]), m1, m2)
                (notBound_BHT, a_BHT, e_BHT, E_BHT) = orbitalElements(x_BHT[N_t_BHT-1], m1, m2)
                #print('a_imp = ', a_imp)
                #print('a_thr = ', a_thr)
                #Energy and a for N body
                E_thr_t = np.zeros(N_t_thr)
                a_thr_t = np.zeros(N_t_thr)
                for i in range(N_t_thr):
                        temp, a_thr_t[i], temp2, E_thr_t[i] = orbitalElements(np.array([x[i,0], x[i,1], x[i,3], x[i,4]]), m1, m2)
                #Energy and a for impulse
                E_imp_t = np.zeros(N_t_imp)
                a_imp_t = np.zeros(N_t_imp)
                for i in range(N_t_imp):
                        temp, a_imp_t[i], temp2, E_imp_t[i] = orbitalElements(x_imp[i], m1, m2)
                #Energy and a for Bahcall et al.
                E_BHT_t = np.zeros(N_t_BHT)
                a_BHT_t = np.zeros(N_t_BHT)
                for i in range(N_t_BHT):
                        temp, a_BHT_t[i], temp2, E_BHT_t[i] = orbitalElements(x_BHT[i], m1, m2)


                
                #PLOTS
                #Plot energy against time
                plt.plot(t_thr*1000.0, E_thr_t*mass_scale()*(length_scale()/time_scale())**2.0, label='N-body')
                plt.plot(t_imp*1000.0, E_imp_t*mass_scale()*(length_scale()/time_scale())**2.0, label='Hyperbolic')
                plt.plot(t_BHT*1000.0, E_BHT_t*mass_scale()*(length_scale()/time_scale())**2.0, label='Double Bahcall et al.')
                #plt.plot(t_thr*1000.0, 0.5*(M_p*mass_scale())*(np.linalg.norm(x[:,5], axis=1)*length_scale()/time_scale())**2.0-0.5*(M_p*mass_scale())*(V_0*length_scale()/time_scale())**2.0, label='Change in kinetic energy of perturber')
                plt.title('Total energy of binary during encounter for initial semimajor axis $10^{}$au and impact parameter $10^{}$au'.format(np.log10(a*length_scale()/au).astype(int), np.log10(b*length_scale()/au).astype(int)), wrap=True)
                plt.xlabel('Time, Myr')
                plt.ylabel('Total energy of binary, J')
                plt.legend()
                plt.show()
                '''
                #Plot perturber speed against time
                plt.plot(t_thr, np.linalg.norm(x[:,5], axis=1))
                plt.title('Perturber speed')
                plt.xlabel('Time')
                plt.ylabel('Perturber speed')
                plt.show()  
                
                #Plot paths of stars
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot(x[:,0,0], x[:,0,1], x[:,0,2])
                ax.plot(x[:,1,0], x[:,1,1], x[:,1,2])
                ax.plot(x[:,2,0], x[:,2,1], x[:,2,2])
                plt.show()
                
                #Plot relative velocity of stars against time
                plt.title('Relative star velocity in x direction')
                plt.plot(t_thr*1000.0, (x[:,3,0]-x[:,4,0])*length_scale()/time_scale(), label='N-body')
                plt.plot(t_imp*1000.0, (x_imp[:,2,0]-x_imp[:,3,0])*length_scale()/time_scale(), label='Impulse')
                plt.xlabel('Time, Myr')
                plt.ylabel('Velocity in x direction, m/s')
                plt.legend()
                plt.show()
                
                plt.title('Relative star velocity in y direction')
                plt.plot(t_thr*1000.0, (x[:,3,1]-x[:,4,1])*length_scale()/time_scale(), label='N-body')
                plt.plot(t_imp*1000.0, (x_imp[:,2,1]-x_imp[:,3,1])*length_scale()/time_scale(), label='Impulse')
                plt.xlabel('Time, Myr')
                plt.ylabel('Velocity in y direction, m/s')
                plt.legend()
                plt.show()
                
                plt.title('Relative star velocity in z direction')
                plt.plot(t_thr*1000.0, (x[:,3,2]-x[:,4,2])*length_scale()/time_scale(), label='N-body')
                plt.plot(t_imp*1000.0, (x_imp[:,2,2]-x_imp[:,3,2])*length_scale()/time_scale(), label='Impulse')
                plt.xlabel('Time, Myr')
                plt.ylabel('Velocity in z direction, m/s')
                plt.legend()
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
                
                                
                #Semimajor axis difference
                a_diff = a_imp - a_thr
                a_frac = a_diff/a_thr
                #Eccentricity difference
                e_diff = e_imp - e_thr
                #print('a_frac = ', a_frac)
                #print('e_diff = ', e_diff)
                
                #print('E_imp =', E_imp)
                #print('E_thr =', E_thr)
                
                #E_frac = (E_imp - E_thr)/E_thr               
                
                print('E_ini =', E_ini)
                delta_E_imp = E_imp - E_ini
                delta_E_thr = E_thr - E_ini
                delta_E_BHT = E_BHT - E_ini
                E_frac = (delta_E_imp - delta_E_thr)/delta_E_thr
                print('delta_E_imp =', delta_E_imp)
                print('delta_E_thr =', delta_E_thr)
                print('delta_E_BHT =', delta_E_BHT)
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

                        
        
        
        
        
        
        
        
        
