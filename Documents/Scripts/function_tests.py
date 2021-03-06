#A script to test some of my functions
import os
os.system("python setup.py build_ext --inplace")

import numpy as np
import matplotlib.pyplot as plt

from orbital_elements import semimajorAxis, orbitalElements
from random_binary import setupRandomBinary, setupRandomBinaryBHT
from evolve_binary import integrateBinary
from scipy.constants import au, parsec

from scipy.constants import G

#Firstly test semimajorAxis and orbitalElements
def orbitalElementsTest():
        #Set up a binary with known orbital elements:
        a = 10**5.0 * au
        e = 0.7
        m1 = 2.0*10.0**30.0
        m2 = 2.0*10.0**30.0

        X = setupRandomBinaryBHT(a, e, m1, m2)
        
        (notBound, a_found, e_found) = orbitalElements(X, m1, m2)
        
        print('Error in a = ', (a_found-a))
        print('Error in e = ', (e_found-e))

        
orbitalElementsTest()     

#Test integrateBinary
def integrateBinaryTest():
        N_t = 10000
        a = 1.0 * au
        e = 0.7
        m1 = 2.0*10.0**30.0
        m2 = 2.0*10.0**30.0
        
        X = np.zeros((N_t,4,3), dtype=float)
        X[0] = setupRandomBinary(a, e, m1, m2)
        t = np.zeros(N_t, dtype=float)
        for i in range(1, N_t):
                #Evolve orbit
                (X[i], dt) = integrateBinary(2, X[i-1], np.array([m1,m2]), n=20)
                #Add time step to time array
                t[i] = t[i-1]+dt
        
        #Total energy
        E = 0.5*m1*np.linalg.norm(X[:,2], axis=1)**2.0 + 0.5*m2*np.linalg.norm(X[:,3], axis=1)**2.0 - G*m1*m2/np.linalg.norm(X[:,0]-X[:,1], axis=1)
                
        plt.plot((X[:,0,0]-X[:,1,0]), (X[:,0,1]-X[:,1,1]))
        plt.show()

        plt.plot(t, (E-E[0])/E[0])
        plt.show()
        
#integrateBinaryTest()