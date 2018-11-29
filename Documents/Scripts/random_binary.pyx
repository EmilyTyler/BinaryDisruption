import random
import numpy as np
cimport numpy as np

from eccentric_anomaly import findEccentricAnomaly

#from scipy.constants import G
from internal_units import *
G = G()

#Function to initialise randomly oriented binary
def setupRandomBinary(double a, double e, double m1, double m2):
        
        cdef double M, E, f, r, n
        cdef np.ndarray x1 = np.zeros(3, dtype=float)
        cdef np.ndarray v1 = np.zeros(3, dtype=float)
        cdef np.ndarray x2 = np.zeros(3, dtype=float)
        cdef np.ndarray v2 = np.zeros(3, dtype=float)
        cdef np.ndarray R = np.zeros(3, dtype=float)
        cdef np.ndarray V = np.zeros(3, dtype=float)

        
        #Randomise mean anomaly
        M = np.random.uniform(0.0, 2.0*np.pi)
        #Find eccentric anomaly from Kepler's equation
        E = findEccentricAnomaly(e, M)
        #Find true anomaly
        f = 2.0*np.arctan(((1.0+e)/(1.0-e))**0.5 * np.tan(0.5*E))
        
        #Randomise true anomaly
        #f = np.random.uniform(0.0, 2.0*np.pi)
        
        #Initial separation
        r = a*(1.0 - e**2.0)/(1.0 + e*np.cos(f))
        
        #Mean motion
        n = np.sqrt(G*(m1+m2)/(a**3.0))
        #Initial coordinates of first star (cartesian)
        x1 = np.array([0.0, 0.0, 0.0])
        #Initial velocity of first star
        v1 = np.array([0.0, 0.0, 0.0])
        #Initial coordinates of second star (cartesian)
        x2 = np.array([r*np.cos(f), r*np.sin(f), 0.0])
        #Initial velocity of second star
        v2 = np.array([- n*a/(np.sqrt(1.0-e**2.0)) * np.sin(f), n*a/(np.sqrt(1.0-e**2.0)) * (e + np.cos(f)), 0.0])
        
        
        #Centre of mass position vector
        R = (m1*x1 + m2*x2)/(m1 + m2)
        #Centre of mass velocity vector
        V = (m1*v1 + m2*v2)/(m1 + m2)
        #Move into centre of mass rest frame
        #x1 -= R
        #x2 -= R
        v1 -= V
        v2 -= V
        
        
        return np.array([x1, x2, v1, v2])

#Function to initialise randomly oriented binary
def setupRandomBinaryBHT(double a, double e, double m1, double m2):
        
        cdef double f, r
        cdef np.ndarray x1 = np.zeros(3, dtype=float)
        cdef np.ndarray v1 = np.zeros(3, dtype=float)
        cdef np.ndarray x2 = np.zeros(3, dtype=float)
        cdef np.ndarray v2 = np.zeros(3, dtype=float)
        cdef np.ndarray R = np.zeros(3, dtype=float)
        cdef np.ndarray V = np.zeros(3, dtype=float)

        e = 0.7
        
        #Randomise mean anomaly
        M = np.random.uniform(0.0, 2.0*np.pi)
        #Find eccentric anomaly from Kepler's equation
        E = findEccentricAnomaly(e, M)
        #Find true anomaly
        f = 2.0*np.arctan(((1.0+e)/(1.0-e))**0.5 * np.tan(0.5*E))
        
        #Initial separation
        #r = a*(1.0 - e**2.0)/(1.0 + e*np.cos(f))
        #Time averaged separation
        #r = a * (1.0 + 0.5*e**2.0)
        #Maximum separation
        #r = a * (1.0 + e)
        #Minimum separation
        r = a * (1.0 - e)   

        #Mean motion
        n = np.sqrt(G*(m1+m2)/(a**3.0))
        #Initial coordinates of first star (cartesian)
        x1 = np.array([0.0, 0.0, 0.0])
        #Initial velocity of first star
        v1 = np.array([0.0, 0.0, 0.0])
        #Initial coordinates of second star (cartesian)
        x2 = np.array([r*np.cos(f), r*np.sin(f), 0.0])
        #Initial velocity of second star
        v2 = np.array([- n*a/(np.sqrt(1.0-e**2.0)) * np.sin(f), n*a/(np.sqrt(1.0-e**2.0)) * (e + np.cos(f)), 0.0])
        #v2 = np.array([0.0, 0.0, 0.0])

        #Centre of mass position vector
        R = (m1*x1 + m2*x2)/(m1 + m2)
        #Centre of mass velocity vector
        V = (m1*v1 + m2*v2)/(m1 + m2)
        #Move into centre of mass rest frame
        x1 -= R
        x2 -= R
        v1 -= V
        v2 -= V
        
        
        return np.array([x1, x2, v1, v2])

