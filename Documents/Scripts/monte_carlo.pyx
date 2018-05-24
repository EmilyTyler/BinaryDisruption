#Draw random numbers from a distribution using a Monte Carlo method

import cython
import random
import numpy as np
cimport numpy as np
from cpython cimport bool

from encounters import calc_b_max, encounter, encounterRate
from scipy.stats import maxwell
from orbital_elements import orbitalElements
from random_binary import setupRandomBinary
from scipy.constants import G

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
def MCEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0):
        #Minimum impact parameter
        cdef double b_min = 0.0
        #Minimum velocity
        cdef double v_min = 10.0**(-2.0)*v_rms
        #Maximum velocity
        cdef double v_max = 10.0**2.0*v_rms
        #Implement encounters
        cdef bool notBound = False
        cdef int N_enc, i, k
        cdef int N_broken = 0
        cdef double b_max
        cdef np.ndarray b, v
        cdef np.ndarray a = np.array([a_0[i] for i in range(N_bin)])
        cdef np.ndarray e = np.array([e_0[i] for i in range(N_bin)])
        for i in range(N_bin):
                #Maximum impact parameter
                #print('a[i] =', a[i])
                b_max = calc_b_max(M_p, v_rms, a[i], m1, m2, prefactor=prefactor)
                #print('b_max =', b_max)
                #Mean number of encounters in time T 
                N_mean = T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max)
                #print('N_mean =', N_mean)
                #Number of encounters in time T 
                N_enc = np.random.poisson(N_mean)
                #print('N_enc =', N_enc)
                #Impact parameters of encounters
                b = draw_b(b_max, N_enc)
                #Relative velocities of encounters
                v = draw_vmaxwellian(v_rms, v_min, v_max, N_enc)
                #print('Implementing',N_enc,'encounters')
                for j in range(N_enc):
                        (notBound, a[i], e[i]) = encounter(m1, m2, v[j], b[j], a[i], e[i], M_p)
                        if notBound:
                                N_broken += 1
                                #print('Binary broken!')
                                a[i] = -1.0
                                e[i] = -1.0
                                break
        return (a, e, N_broken)

#Jiang and Tremaine simulation of encounters of N_bin binaries over time T
def JTEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a_0, np.ndarray[double, ndim=1] e_0, int N_bin, double prefactor = 1.0):
        cdef double rho = M_p * n_p
        cdef np.ndarray m = np.array([m1, m2])
        #Implement encounters
        cdef bool notBound = False
        cdef int N_enc, i, k
        cdef int N_broken = 0
        cdef double b_max
        cdef np.ndarray b, v, X
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
                if notBound:
                        a[i] = -1.0
                        e[i] = -1.0
                        N_broken += 1                                
        return (a, e, N_broken)



        
