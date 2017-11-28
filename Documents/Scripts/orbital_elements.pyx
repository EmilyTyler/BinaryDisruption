#Contains functions to find the orbital elements of a binary star given its position and velocity vectors

import numpy as np
cimport numpy as np

from scipy.constants import G

#To find the semi-major axis of a binary star
def semimajorAxis(np.ndarray[double, ndim=2] X, double m1, double m2):
        cdef double R, V, E, a
        R = np.linalg.norm(X[0]-X[1])
        V = np.linalg.norm(X[2]-X[3])
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G/R)
        #Semi-major axis
        a = G*m1*m2/(2.0*abs(E))
        return a



def orbitalElements(np.ndarray[double, ndim=2] X, double m1, double m2):
        cdef np.ndarray x = np.zeros(3, dtype=float)
        cdef np.ndarray v = np.zeros(3, dtype=float)
        cdef np.ndarray L = np.zeros(3, dtype=float)
        cdef double R, V, E, a, e
        
        x = X[0] - X[1]
        v = X[2] - X[3]
        R = np.linalg.norm(x)
        V = np.linalg.norm(v)
        
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G/R)       
        #Total angular momentum
        L = m1*m2/(m1+m2)*np.cross(x,v)       
        #Semi-major axis
        a = G*m1*m2/(2.0*abs(E))      
        #Eccentricity
        e = np.sqrt(1.0 + 2.0*(m1+m2)*np.dot(L,L)*E/(G**2.0*(m1*m2)**3.0))
        
        return((E >= 0.0), a, e)
       



        
        
