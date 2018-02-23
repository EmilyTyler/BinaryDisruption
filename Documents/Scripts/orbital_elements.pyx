#Contains functions to find the orbital elements of a binary star given its position and velocity vectors

import numpy as np
cimport numpy as np

from scipy.constants import G

#To find the semi-major axis of a binary star
def semimajorAxis(np.ndarray[double, ndim=2] X, double m1, double m2):
        cdef np.ndarray x = np.zeros(3, dtype=float)
        cdef np.ndarray v = np.zeros(3, dtype=float)
        cdef double R, V, E, a
        R = np.sqrt(x[0]**2.0 + x[1]**2.0 + x[2]**2.0)
        V = np.sqrt(v[0]**2.0 + v[1]**2.0 + v[2]**2.0)
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G/R)
        #Semi-major axis
        a = G*m1*m2/(2.0*abs(E))
        return a



def orbitalElements(np.ndarray[double, ndim=2] X, double m1, double m2):
        cdef np.ndarray x = np.zeros(3, dtype=float)
        cdef np.ndarray v = np.zeros(3, dtype=float)
        cdef np.ndarray L = np.zeros(3, dtype=float)
        cdef double R, V, E, a, e, L_squared
        
        x = X[0] - X[1]
        v = X[2] - X[3]
        R = np.sqrt(x[0]**2.0 + x[1]**2.0 + x[2]**2.0)
        V = np.sqrt(v[0]**2.0 + v[1]**2.0 + v[2]**2.0)
        
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G/R)       
        #Total angular momentum
        L = m1*m2/(m1+m2)*np.array([x[1]*v[2]-x[2]*v[1], x[2]*v[0]-x[0]*v[2], x[0]*v[1]-x[1]*v[0]])
        L_squared = L[0]**2.0 + L[1]**2.0 + L[2]**2.0
        #Semi-major axis
        a = G*m1*m2/(2.0*abs(E))      
        #Eccentricity
        e = np.sqrt(1.0 + 2.0*(m1+m2)*L_squared*E/(G**2.0*(m1*m2)**3.0))
        
        return((E >= -10.0**(-4.0)), a, e)
       



        
        
