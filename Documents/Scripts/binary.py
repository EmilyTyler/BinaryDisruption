import numpy as np
import math
import random

from eccentric_anomaly import findEccentricAnomaly


#Initialise variables
#Semi-major axis, pc
a = 0.1
#Eccentricity
e = 0.7

#RMS of Maxwellian velocity distribution, km/s
v_rms = 100.0

#Convert to SI
a = a * 3.086**(16.0)
v_rms = v_rms * 1000.0

#Function to find perturber velocity
def  relativeVelocity():
        #Velocity is in z-direction
        v_rel = [0,0,v_rms]
        return v_rel


#Initialise binary
#Find mean anomaly
M = random.uniform(0.0, 2.0*math.pi)
print M
