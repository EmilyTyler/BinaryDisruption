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
#Randomise mean anomaly
M = random.uniform(0.0, 2.0*math.pi)
#Find eccentric anomaly from Kepler's equation
E = findEccentricAnomaly(e, M)
#Find true anomaly
f = 2.0*math.atan(((1.0+e)/(1.0-e))**0.5 * math.tan(0.5*E))


