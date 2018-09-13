#Testing algorithms to solve Kepler's equation
import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from eccentric_anomaly import findEccentricAnomaly, findEccentricAnomaly2

e = 0.7
M = np.random.uniform(0.0, 2.0*np.pi)

E1, count1 = findEccentricAnomaly(e, M)
E2, count2 = findEccentricAnomaly2(e, M)

print('Method 1:')
print('E =', E1)
print('count =', count1)
print('Method 2:')
print('E =', E2)
print('count =', count2)