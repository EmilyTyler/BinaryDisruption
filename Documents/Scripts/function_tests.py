#A script to test some of my functions

import numpy as np
from orbital_elements import semimajorAxis, orbitalElements

G = 6.67 * 10.0**(-11.0)

#Firstly test semimajorAxis and orbitalElements
def orbitalElementsTest():
        #Set up a binary with known orbital elements:
        a = 0.1 * 3.086*10.0**16.0
        e = 0.7
        f = 0.0
        I = 2.0
        Omega = 1.0
        omega = -1.5
        m1 = 2.0*10.0**30.0
        m2 = 2.0*10.0**30.0
        
        r = a*(1.0 - e**2.0)/(1.0 + e*np.cos(f))
        x1 = np.array([0.0, 0.0, 0.0])
        v1 = np.array([0.0, 0.0, 0.0])
        #Mean motion
        n = np.sqrt(G*(m1+m2)/(a**3.0))
        x2 = np.array([r*np.cos(f), r*np.sin(f), 0.0])
        v2 = np.array([- n*a/(np.sqrt(1.0-e**2.0)) * np.sin(f), n*a/(np.sqrt(1.0-e**2.0)) * (e + np.cos(f)), 0.0])
        X = np.array([x1, x2, v1, v2])
        #Rotate
        R1 = np.array([[np.cos(omega), -np.sin(omega), 0.0],
                       [np.sin(omega), np.cos(omega), 0.0],
                       [0.0, 0.0, 1.0]])
        R2 = np.array([[1.0, 0.0, 0.0],
                       [0.0, np.cos(I), -np.sin(I)],
                       [0.0, np.sin(I), np.cos(I)]])
        R3 = np.array([[np.cos(Omega), -np.sin(Omega), 0.0],
                       [np.sin(Omega), np.cos(Omega), 0.0],
                       [0.0, 0.0, 1.0]])
        R = np.dot(np.dot(R3, R2), R1)
        X[1] = np.transpose(np.dot(R, np.transpose(X[1])))
        X[3] = np.transpose(np.dot(R, np.transpose(X[3])))
        
        (a_found, e_found, f_found, I_found, Omega_found, omega_found) = orbitalElements(X, m1, m2)
        
        print('Error in a = ', (a_found-a))
        print('Error in e = ', (e_found-e))
        print('Error in f = ', (f_found-f))
        print('Error in I = ', (I_found-I))
        print('Error in Omega = ', (Omega_found-Omega))
        print('Error in omega = ', (omega_found-omega))
        
orbitalElementsTest()     