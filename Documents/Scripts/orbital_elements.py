#Contains functions to find the orbital elements of a binary star given its position and velocity vectors

import numpy as np

G = 6.67*10.0**(-11.0)

#To find the semi-major axis of a binary star
def semimajorAxis(X, m1, m2):
        R = np.linalg.norm(X[0]-X[1])
        V = np.linalg.norm(X[2]-X[3])
        a = (2.0/R - V**2.0/(G*(m1+m2)))**(-1.0)
        return a



def orbitalElements(X, m1, m2):
        R = np.linalg.norm(X[0]-X[1])
        V = np.linalg.norm(X[2]-X[3])
        R_dot_Rdot = np.dot((X[0]-X[1]),(X[2]-X[3]))
        h = np.cross((X[0]-X[1]),(X[2]-X[3]))
        if abs(V**2.0 - np.dot(h,h)/(R**2.0)) > 10.0**(-8.0):
                R_dot = np.sign(R_dot_Rdot)*np.sqrt(V**2.0 - np.dot(h,h)/(R**2.0))
        else:
                R_dot = 0.0
                        
        
        #Semi-major axis
        a = (2.0/R - V**2.0/(G*(m1+m2)))**(-1.0)
        #Eccentricity
        e = np.sqrt(1.0 - np.dot(h,h)/(G*(m1+m2)*a))
        #Inclination
        I = np.arccos(h[2]/np.linalg.norm(h))
        #True anomaly
        f = np.arcsin(a*(1.0-e**2.0)*R_dot/(np.linalg.norm(h)*e))
        #
        Omega = np.arctan(-h[0]/h[1])
        print(h[1])
        print(-h[0]/h[1])
        print(Omega)
        if I != 0.0:
                omega_plus_f = np.arcsin((X[0,2]-X[1,2])/R*np.sin(I))
                omega = omega_plus_f - f
        else:
                omega = 0.0
                
        print(omega)
                
        return(a, e, f, I, Omega, omega)