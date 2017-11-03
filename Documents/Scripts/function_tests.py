#A script to test some of my functions

import numpy as np
from orbital_elements import semimajorAxis, orbitalElements
from random_binary import setupRandomBinary

G = 6.67 * 10.0**(-11.0)

#Firstly test semimajorAxis and orbitalElements
def orbitalElementsTest():
        #Set up a binary with known orbital elements:
        a = 0.1 * 3.086*10.0**16.0
        e = 0.7
        m1 = 2.0*10.0**30.0
        m2 = 2.0*10.0**30.0

        X = setupRandomBinary(a, e, m1, m2)
        
        (a_found, e_found) = orbitalElements(X, m1, m2)
        
        print('Error in a = ', (a_found-a))
        print('Error in e = ', (e_found-e))

        
orbitalElementsTest()     
