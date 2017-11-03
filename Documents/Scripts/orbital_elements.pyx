#Contains functions to find the orbital elements of a binary star given its position and velocity vectors

import numpy as np

G = 6.67*10.0**(-11.0)

#To find the semi-major axis of a binary star
def semimajorAxis(X, m1, m2):
        R = np.linalg.norm(X[0]-X[1])
        V = np.linalg.norm(X[2]-X[3])
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G/R)
        #Semi-major axis
        a = G*m1*m2/(2.0*abs(E))
        return a



def orbitalElements(X, m1, m2):
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
                               
        return(a, e)
        
        
