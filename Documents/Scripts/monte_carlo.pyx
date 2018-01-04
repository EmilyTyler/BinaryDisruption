#Draw random numbers from a distribution using a Monte Carlo method


import random
import numpy as np
cimport numpy as np
from cpython cimport bool

from encounters import calc_b_max
from encounters import encounterRate
from encounters import encounter

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

#Draw impact parameter from a distribution linear in b
def draw_b(double b_max, int N):
        return b_max*np.sqrt(np.random.uniform(0.0, 1.0, size=N))

#Monte Carlo simulation of encounters of N_bin binaries over time T
def MCEncounters(double v_rms, double n_p, double T, double m1, double m2, double M_p, np.ndarray[double, ndim=1] a, np.ndarray[double, ndim=1] e, int N_bin):
        #Minimum impact parameter
        cdef double b_min = 10.0**(-2.0)*(np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
        #Minimum velocity
        cdef double v_min = 10.0**(-2.0)*v_rms
        #Maximum velocity
        cdef double v_max = 10.0**2.0*v_rms
        #Implement encounters
        cdef bool notBound = False
        cdef int N_enc
        cdef double a_new, e_new, b_max
        cdef np.ndarray b, v
        for i in range(N_bin):
                #Maximum impact parameter
                b_max = calc_b_max(M_p, v_rms, a[i], m1, m2)
                #Total number of encounters in time T 
                N_enc = int(T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max))
                #Impact parameters of encounters
                b = draw_b(b_max, N_enc)
                #Relative velocities of encounters
                v = draw_maxwellian(v_rms, v_min, v_max, N_enc)
                for j in range(N_enc):
                        (notBound, a_new, e_new) = encounter(m1, m2, v[j], b[j], a[i], e[i], M_p)
                        if notBound:
                                print('Binary broken!')
                                a[i] = 0.0
                                e[i] = 0.0
                                break
                        a[i] = a_new
                        e[i] = e_new
        return (a, e)



        
