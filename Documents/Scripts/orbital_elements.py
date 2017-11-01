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
        x = X[1] - X[0]
        v = X[3] - X[2]
        R = np.linalg.norm(x)
        V = np.linalg.norm(v)
        R_dot_Rdot = np.dot(x,v)
        h = np.cross(x,v)
        if (-10.0**(-10.0)< (V**2.0 - np.dot(h,h)/(R**2.0)) < 0.0):
                R_dot = 0.0
        else:
                R_dot = np.sign(R_dot_Rdot)*np.sqrt(V**2.0 - np.dot(h,h)/(R**2.0))            
        #Semi-major axis
        a = (2.0/R - V**2.0/(G*(m1+m2)))**(-1.0)
        #Eccentricity
        e = np.sqrt(1.0 - np.dot(h,h)/(G*(m1+m2)*a))
        #Inclination
        I = np.arccos(h[2]/np.linalg.norm(h))
        #True anomaly
        f = np.arcsin(a*(1.0-e**2.0)*R_dot/(np.linalg.norm(h)*e))
        if abs(np.cos(f) - (a*(1.0-e**2.0)/R - 1.0)/e) > 10.0**(-10.0):
                f = np.sign(f)*np.pi - f
                if abs(np.cos(f) - (a*(1.0-e**2.0)/R - 1.0)/e) > 10.0**(-10.0):
                        print('f not found')
        #Longitude of ascending node
        if I != 0.0:            
                Omega = np.arcsin(np.sign(h[2])*h[0]/(np.linalg.norm(h)*np.sin(I)))
                if abs(np.cos(Omega) + np.sign(h[2])*h[1]/(np.linalg.norm(h)*np.sin(I))) > 10.0**(-10.0):
                        Omega = np.sign(Omega)*np.pi - Omega
                        if abs(np.cos(Omega) + np.sign(h[2])*h[1]/(np.linalg.norm(h)*np.sin(I))) > 10.0**(-10.0):
                                print('Omega not found')            
        else:
                Omega = 0.0
        #Argument of pericentre
        omega_plus_f = np.arcsin(x[2]/R*np.sin(I))
        print(abs(np.cos(omega_plus_f) - ((x[0])/R + np.sin(Omega)*np.sin(omega_plus_f)*np.cos(I))/np.cos(Omega)))
        if abs(np.cos(omega_plus_f) - ((x[0])/R + np.sin(Omega)*np.sin(omega_plus_f)*np.cos(I))/np.cos(Omega)) > 10.0**(-10.0):
                print('here')
                omega_plus_f = np.sign(omega_plus_f)*np.pi - omega_plus_f
                print(abs(np.cos(omega_plus_f) - ((x[0])/R + np.sin(Omega)*np.sin(omega_plus_f)*np.cos(I))/np.cos(Omega)))
                print(abs(np.sin(omega_plus_f) - x[2]/R*np.sin(I)))
                if abs(np.cos(omega_plus_f) - ((x[0])/R + np.sin(Omega)*np.sin(omega_plus_f)*np.cos(I))/np.cos(Omega)) > 10.0**(-10.0):
                        print('omega_plus_f not found')
        print(omega_plus_f - f)                
        omega = (omega_plus_f - f) % (2.0*np.pi)
        if omega > np.pi:
                omega = omega - 2.0*np.pi
                               
        return(a, e, f, I, Omega, omega)
        
        
