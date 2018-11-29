#Draw random numbers from a distribution using a Monte Carlo method

import cython
import random
import numpy as np
cimport numpy as np
from cpython cimport bool

from encounters import calc_b_max, encounter, encounterRate, impulseEncounter, BHTEncounter
from scipy.stats import maxwell
from orbital_elements import orbitalElements
from random_binary import setupRandomBinary
from scipy.constants import parsec, au, mega, year
from internal_units import *
G = G()

#Draw velocities from a Maxwellian distribution
#Uses rejection method from Numerical Recipes, Press et al. section 7.3.6
def draw_maxwellian(double v_rms, double v_min, double v_max, int N):      
        #Accepted points
        cdef np.ndarray accepted = np.array([])
        #Area 
        cdef double area = (v_max-v_min)*4.0*np.exp(-1.0)*(2.0*np.pi)**(-0.5)/v_rms
        
        cdef double u, x, y
        while accepted.size < N:                
                u = random.uniform(0.0, area)
                x = maxwellianX_from_area(u, v_rms, v_min)
                y = random.uniform(0.0, maxwellianComparison(x, v_rms, v_min, v_max))
                if y < maxwellianPdf(x, v_rms):
                        accepted = np.append(accepted, x)                      
        return accepted


def maxwellianPdf(double x, double v_rms):     
        return 4.0*np.pi*x**2.0*np.exp(-x**2.0/(2.0*v_rms**2.0))*((2.0*np.pi*v_rms**2.0))**(-3.0/2.0)

def maxwellianComparison(double x, double v_rms, double v_min, double v_max):      
        cdef double answer = 0.0
        if v_min<x<v_max:
                answer = 4.0*np.exp(-1.0)*(2.0*np.pi)**(-0.5)/v_rms
        return answer

def maxwellianX_from_area(double A, double v_rms, double v_min):      
        return (2.0*np.pi)**0.5*v_rms*np.exp(1.0)*A/4.0 + v_min

#Draw velocities from a Maxwellian times v distribution
def draw_vmaxwellian(double v_rms, double v_min, double v_max, int N):
        #Accepted points
        cdef np.ndarray accepted = np.array([])
        #Total area under comparison function
        cdef double area = (v_max - v_min)*3.0**(3.0/2.0)/(2.0*v_rms)*np.exp(-3.0/2.0)
        
        cdef double u, x, y
        while accepted.size < N:
                u = random.uniform(0.0, area)
                x = vmaxwellianX_from_area(u, v_rms, v_min)
                y = random.uniform(0.0, vmaxwellianComparison(x, v_rms, v_min, v_max))
                if y < vmaxwellianPdf(x, v_rms):
                        accepted = np.append(accepted, x)
        return accepted

def vmaxwellianPdf(double x, double v_rms):
        return x**3.0/(2.0*v_rms**4.0)*np.exp(-x**2.0/(2.0*v_rms**2.0))

def vmaxwellianComparison(double x, double v_rms, double v_min, double v_max):
        cdef double answer = 0.0
        if v_min<x<v_max:
                answer = 3.0**(3.0/2.0)/(2.0*v_rms)*np.exp(-3.0/2.0)
        return answer

def vmaxwellianX_from_area(double A, double v_rms, double v_min):
        return v_min + 2.0*v_rms*A*np.exp(3.0/2.0)/3.0**(3.0/2.0)


#Draw impact parameter from a distribution linear in b
def draw_b(double b_max, int N):
        return b_max*np.sqrt(np.random.uniform(0.0, 1.0, size=N))

#Monte Carlo simulation of encounters of N_bin binaries over time T
def MCEncounters_new(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0, double a_T = parsec):
        #Minimum impact parameter
        cdef double b_min = 0.0
        #Maximum maximum impact parameter
        cdef double b_max_max = calc_b_max(M_p, v_rms, a_T, m1, m2, prefactor=prefactor)
        print('b_max_max, pc =', b_max_max*length_scale()/parsec)
        #Minimum velocity
        cdef double v_min = 10.0**(-2.0)*v_rms
        #Maximum velocity
        cdef double v_max = 10.0**2.0*v_rms
        #Mean number of encounters in time T 
        N_mean = T*encounterRate(n_p, v_rms, b_min, b_max_max, v_min, v_max)
        #print('N_mean =', N_mean)
        #Number of encounters
        cdef np.ndarray N_enc = np.random.poisson(N_mean, size=N_bin)
        #cdef np.ndarray N_enc_actual = np.zeros(N_bin, dtype=int)
        #Implement encounters
        cdef bool notBound = False
        cdef int i, j
        cdef int N_broken = 0
        cdef double b_max, v
        cdef np.ndarray b
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        for i in range(N_bin):
                #Impact parameters of encounters
                b = draw_b(b_max_max, N_enc[i])
                #Maximum impact parameter
                b_max = calc_b_max(M_p, v_rms, a[i], m1, m2, prefactor=prefactor)
                #Implement encounters
                for j in range(N_enc[i]):
                        if (b[j]<=b_max):
                                #Draw velocity
                                v = draw_vmaxwellian(v_rms, v_min, v_max, 1)[0]
                                #Encounter
                                (notBound, a[i], e[i]) = impulseEncounter(m1, m2, v, b[j], a[i], e[i], M_p)
                                #N_enc_actual[i] += 1
                                if (notBound or a[i]>=a_T):
                                        N_broken += 1
                                        #print('Binary broken!')
                                        a[i] = -1.0
                                        e[i] = -1.0
                                        break
                                #Update maximum impact parameter
                                b_max = calc_b_max(M_p, v_rms, a[i], m1, m2, prefactor=prefactor)
                #if (a[i]>0.0):
                        #print('N_enc =', N_enc_actual[i])
        return (a, e, N_broken)

#Monte Carlo simulation of encounters of N_bin binaries over time T
def MCEncounters_t(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0):
        #Minimum impact parameter
        cdef double b_min = 0.0
        #Minimum velocity
        cdef double v_min = 10.0**(-2.0)*v_rms
        #Maximum velocity
        cdef double v_max = 10.0**2.0*v_rms
        #Implement encounters
        cdef bool notBound = False
        cdef int i, j
        cdef int N_broken = 0
        cdef double b_max, v, t, E_fin
        cdef np.ndarray b
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        cdef double a_MW = 1000*parsec/length_scale()
        for i in range(N_bin):
                #Time passed
                t = 0.0
                #Implement encounters
                while t <= T:
                        #Maximum impact parameter
                        b_max = calc_b_max(M_p, v_rms, a[i], m1, m2, prefactor=prefactor)
                        #print('b_max, pc =', b_max*length_scale()/parsec)
                        #Encounter rate 
                        rate = encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max)
                        #Increment time passed
                        t += np.random.exponential(1.0/rate)
                        #print('t/T =', t/T)
                        #Draw velocity
                        v = draw_vmaxwellian(v_rms, v_min, v_max, 1)[0]
                        #Draw impact parameter
                        b = draw_b(b_max, 1)
                        #Encounter
                        #print('a[i], pc =', a[i]*length_scale()/parsec)
                        (notBound, a[i], e[i], E_fin) = impulseEncounter(m1, m2, v, b, a[i], e[i], M_p)
                        #print('E_fin =', E_fin)
                        #print('notBound =', notBound)
                        if (notBound or a[i]>a_MW):
                                N_broken += 1
                                #print('Binary broken!')
                                a[i] = -1.0
                                e[i] = -1.0
                                break
                        
        return (a, e, N_broken)

#Jiang and Tremaine simulation of encounters of N_bin binaries over time T
def JTEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0, double a_T=10.0**10.0*au):
        cdef double rho = M_p * n_p
        cdef np.ndarray m = np.array([m1, m2])
        #Implement encounters
        cdef bool notBound = False
        cdef int i, j, k
        cdef int N_broken = 0
        cdef double b_max
        cdef np.ndarray X
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        for i in range(N_bin):
                #Maximum impact parameter
                b_max = a_0[i]/2.0
                #Setup binary
                X = setupRandomBinary(a_0[i], e_0[i], m1, m2)
                L = b_max*v_rms**2.0/(G*(m+M_p))
                #First order diffusion coefficient for parallel velocity
                D_v_para_1 = 4.0*(2.0*np.pi)**0.5*G**2.0*(m+M_p)*rho*np.log(L)*(X[2:,0]**2.0+X[2:,1]**2.0+X[2:,2]**2.0)**0.5/(3.0*v_rms**3.0)
                #Second order diffusion coefficient for perpendicular velocity
                D_v_perp_2 = 16.0*(2.0*np.pi)**0.5*G**2.0*M_p*rho*np.log(L)/(3.0*v_rms)
                #First order diffusion coefficient
                D_v_1 = np.zeros((2,3), dtype=float)
                for j in range(2):
                        #Calculate diffusion coefficients
                        v_j = (X[j+2,0]**2.0+X[j+2,1]**2.0+X[j+2,2]**2.0)**0.5
                        for k in range(3):
                                D_v_1[j,k] = X[j+2,k]/v_j*D_v_para_1[j]
                #Calculate mean and standard deviations
                mean = D_v_1 * T
                std_dev = np.sqrt(D_v_perp_2*T/2.0)
                #Draw velocity changes
                dv = np.zeros((2,3))
                for j in range(2):
                        for k in range(3):
                                dv[j,k] = np.random.normal(mean[j,k], std_dev[j])
                #Add velocity changes to velocities
                X[2:] += dv
                #Close binary
                notBound, a[i], e[i] = orbitalElements(X, m1, m2)
                if notBound or a[i]>=a_T:
                        a[i] = -1.0
                        e[i] = -1.0
                        N_broken += 1                                
        return (a, e, N_broken)

#Weinberg et al encounters
def WSWEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double a_T=10.0**10.0*au):
        cdef int i, j, k, N_t, N_enc
        cdef double dt, N_mean, delta_E, a_new, e_new, q, E_new
        cdef np.ndarray t
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        cdef double M_b = m1 + m2
        cdef np.ndarray E = -G*M_b/(2.0*a)
        #Check we're in the single kick regime
        #print('This should be much greater than 1: ', np.min((M_b/M_p)**2.0*0.1*a*v_rms**2.0/(G*M_b)))
        #Period of binaries
        cdef np.ndarray P = 2.0*np.pi*(a**3.0/(G*M_b))**0.5
        #Maximum impact parameters
        cdef np.ndarray b_max = v_rms*P/(2.0*np.pi)
        #I HAVE TO USE MY OWN B_MAX, THIS IS FAR TOO LARGE
        #cdef np.ndarray b_max = np.array([calc_b_max(M_p, v_rms, a[i], m1, m2, delta=10.0**(-6.0)) for i in range(N_bin)])
        #print('b_max/pc =', b_max/parsec)
        #Impact parameter where encounters at b>b_FP are treated with FP
        cdef np.ndarray b_FP = (16.0*G*M_b/(3.0*0.1*a*v_rms**2.0))**0.5*(M_p/M_b)*a
        #First order diffusion coefficients
        cdef np.ndarray epsilon_1 = 8.0*(2.0*np.pi)**0.5*n_p*P*G**2.0*M_p**2.0*np.log(b_max/b_FP)/v_rms
        #Second order diffusion coefficients
        cdef np.ndarray epsilon_2 = np.sqrt(2.0*epsilon_1*G*M_b/(3.0*a))      
        cdef int N_broken = 0
        for i in range(N_bin):
                #Rate of encounters
                rate = encounterRate(n_p, v_rms, 0.0, b_max[i], 0.01*v_rms, 100.0*v_rms)
                #Timestep
                dt = np.min([np.max([(-0.03*E[i])**2.0*P[i]/epsilon_2[i]**2.0, P[i], 100.0*rate**(-1.0)]), T])
                #
                #print('dt*rate =', dt*rate)
                #Number of timesteps
                #print('T =', T)
                #print('dt =', dt)
                N_t = int(T/dt)
                #Adjust timestep
                dt = T/N_t
                #Time array
                t = np.array([l*dt for l in range(N_t)])
                for j in range(N_t):
                        #Reset energy change
                        delta_E = 0.0
                        #Number of encounters
                        N_enc = np.random.poisson(rate*dt)
                        #Impact parameters
                        bs = draw_b(b_max[i], N_enc)
                        #Impact parameters in the catastrophic regime
                        b_cat = bs[np.where(bs<b_FP[i])]
                        #Draw velocities
                        v = draw_vmaxwellian(v_rms, 0.01*v_rms, 100.0*v_rms, np.size(b_cat))
                        for k in range(np.size(b_cat)):
                                #print('Catastrophic encounter!')                               
                                notBound, a_new, e_new = impulseEncounter(m1, m2, v[k], b_cat[k], a[i], e[i], M_p)   
                                #Energy change
                                delta_E += -G*M_b/(2.0*a_new) + G*M_b/(2.0*a[i])        
                        #Implement FP encounters
                        #Energy change
                        delta_E += epsilon_1[i]*dt/P[i] + np.random.normal()*epsilon_2[i]*(dt/P[i])**0.5
                        #New energy:
                        E_new = delta_E - G*M_b/(2.0*a[i])
                        #Check it's bound
                        if E_new>=0.0 or a_new>=a_T:
                                a[i] = -1.0
                                e[i] = -1.0
                                N_broken += 1
                                break
                        else:
                                #Set new semi-major axis
                                a[i] = -G*M_b/(2.0*E_new)
        return a, e, N_broken
                        
#Monte Carlo simulation of encounters of N_bin binaries over time T
def ClosestEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0, double a_T = parsec, int N_closest=100):
        #Minimum impact parameter
        cdef double b_min = 0.0
        #Maximum maximum impact parameter
        cdef double b_max_max = calc_b_max(M_p, v_rms, a_T, m1, m2, prefactor=prefactor)
        #Minimum velocity
        cdef double v_min = 10.0**(-2.0)*v_rms
        #Maximum velocity
        cdef double v_max = 10.0**2.0*v_rms
        #Mean number of encounters in time T 
        N_mean = T*encounterRate(n_p, v_rms, b_min, b_max_max, v_min, v_max)
        #Number of encounters
        cdef np.ndarray N_enc = np.random.poisson(N_mean, size=N_bin)
        #Implement encounters
        cdef bool notBound = False
        cdef int i, j
        cdef int N_broken = 0
        cdef np.ndarray b, v
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        for i in range(N_bin):
                #Impact parameters of encounters
                b = draw_b(b_max_max, N_enc[i])
                #Find N_closest closest encounters
                b = b[np.where(b < np.sort(b)[N_closest])]
                #Draw velocity
                v = draw_vmaxwellian(v_rms, v_min, v_max, N_closest)
                #Implement encounters
                for j in range(N_closest):
                        #Encounter
                        (notBound, a[i], e[i]) = impulseEncounter(m1, m2, v[j], b[j], a[i], e[i], M_p)
                        if (notBound or a[i]>=a_T):
                                N_broken += 1
                                #print('Binary broken!')
                                a[i] = -1.0
                                e[i] = -1.0
                                break
        return (a, e, N_broken)              

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



        
