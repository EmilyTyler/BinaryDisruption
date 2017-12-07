#Draw random numbers from a distribution using a Monte Carlo method
#Uses rejection method from Numerical Recipes, Press et al. section 7.3.6

import random
import numpy as np
cimport numpy as np


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


def maxwellianPdf(x, v_rms):     
        return 4.0*np.pi*x**2.0*np.exp(-x**2.0/(2.0*v_rms**2.0))*(2.0*np.pi*v_rms**2.0)**(-3.0/2.0)

def maxwellianComparison(x, v_rms, v_min, v_max):      
        if v_min<x<v_max:
                answer = 4.0*np.exp(-1.0)*(2.0*np.pi)**(-0.5)/v_rms
        else:
                answer = 0.0
        return answer

def maxwellianX_from_area(A, v_rms, v_min):      
        return (2.0*np.pi)**0.5*v_rms*np.exp(1.0)*A/4.0 + v_min




        
