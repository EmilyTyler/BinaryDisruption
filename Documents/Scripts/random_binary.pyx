import math
import random
import numpy as np
cimport numpy as np

from eccentric_anomaly import findEccentricAnomaly

G = 6.67 * 10.0**(-11.0)

# cython profile=True

#Function to initialise randomly oriented binary
def setupRandomBinary(double a, double e, double m1, double m2):
        
        cdef double M, E, f, r, n
        cdef np.ndarray x1 = np.zeros(3, dtype=float)
        cdef np.ndarray v1 = np.zeros(3, dtype=float)
        cdef np.ndarray x2 = np.zeros(3, dtype=float)
        cdef np.ndarray v2 = np.zeros(3, dtype=float)
        
        #Randomise mean anomaly
        M = random.uniform(0.0, 2.0*math.pi)
        #Find eccentric anomaly from Kepler's equation
        E = findEccentricAnomaly(e, M)
        #Find true anomaly
        f = 2.0*math.atan(((1.0+e)/(1.0-e))**0.5 * math.tan(0.5*E))
        #Initial separation
        r = a*(1.0 - e**2.0)/(1.0 + e*math.cos(f))
        #Mean motion
        n = math.sqrt(G*(m1+m2)/(a**3.0))
        #Initial coordinates of first star (cartesian)
        x1 = np.array([0.0, 0.0, 0.0])
        #Initial velocity of first star
        v1 = np.array([0.0, 0.0, 0.0])
        #Initial coordinates of second star (cartesian)
        x2 = np.array([r*math.cos(f), r*math.sin(f), 0.0])
        #Initial velocity of second star
        v2 = np.array([- n*a/(math.sqrt(1.0-e**2.0)) * math.sin(f), n*a/(math.sqrt(1.0-e**2.0)) * (e + math.cos(f)), 0.0])
        
        return np.array([x1, x2, v1, v2])