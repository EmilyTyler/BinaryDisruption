#Draw random numbers from a distribution using a Monte Carlo method
#Uses rejection method from Numerical Recipes, Press et al. section 7.3.6

import random
import numpy as np
from matplotlib import pyplot as plt

def rejection(pdf, comparison, x_from_area, area, N):
        
        #Accepted points
        accepted = []
        for i in range(N):
                
                u = random.uniform(0.0, area)
                x = x_from_area(u)
                y = random.uniform(0.0, comparison(x))
                if y < pdf(x):
                        accepted.append(x)
                        
        return accepted


def maxwellianPdf(x):
        
        return 4.0*np.pi*x**2.0*np.exp(-x**2.0/(2.0*v_rms**2.0))*(2.0*np.pi*v_rms**2.0)**(-3.0/2.0)

def maxwellianComparison(x):
        
        if v_min<x<v_max:
                answer = 4.0*np.exp(-1.0)*(2.0*np.pi)**(-0.5)/v_rms
        else:
                answer = 0.0
        return answer

def maxwellianX_from_area(A):
        
        return (2.0*np.pi)**0.5*v_rms*np.exp(1.0)*A/4.0 + v_min

v_rms = 10.0**5.0
v_max = 10.0**6.0
v_min = 10.0**(-3.0)
area = (v_max-v_min)*4.0*np.exp(-1.0)*(2.0*np.pi)**(-0.5)/v_rms
N = 10**7

max_dist = rejection(maxwellianPdf, maxwellianComparison, maxwellianX_from_area, area, N)

#Create velocity bins
N_v = 1000
dlogv = (np.log(v_max) - np.log(v_min))/N_v
v_bins = np.array([v_min*np.exp(dlogv*i) for i in range(N_v)])

#Number of points in each bin
p = np.zeros(N_v, dtype = int)
for x in max_dist:
        i = int(np.floor(np.log(x/v_min)/dlogv))
        p[i] += 1
        
     
plt.plot(v_bins/1000.0, p)
plt.xlabel('v, km/s')
plt.ylabel('pdf')
plt.show()


        