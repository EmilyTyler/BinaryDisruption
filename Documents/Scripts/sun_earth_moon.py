import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import year, au
from evolve_binary import integrateBinary

#To test 3 body integration
#Sun, Earth, Moon

#Time array
t = np.array([0.0])
#Finish time
t_end = 1.0 * year
#Initial positions and velocities
x = np.array([[[0.0, 0.0, 0.0],
              [0.0, 1.0*au, 0.0],
              [0.0, 1.0*au+3.85*10.0**8.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 3.0*10.0**4.0],
              [0.0, 0.0, 3.0*10.0**4.0+1.022*10.0**3.0]]])
#Masses
m = np.array([2.0*10.0**30.0, 5.972*10.0**24.0, 7.348*10.0**22.0])

#Initialise counter
i = 1
while t[i-1] < t_end:
        (x_new, dt) = integrateBinary(3, x[i-1], m, n=20)
        x = np.append(x, [x_new], axis=0)
        t = np.append(t, t[i-1]+dt)
        #Increment counter
        i += 1
        
#Plot paths in Sun's rest frame
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:,0,0], x[:,0,1], x[:,0,2])
ax.plot(x[:,1,0], x[:,1,1], x[:,1,2])
ax.plot(x[:,2,0], x[:,2,1], x[:,2,2])
plt.show()

#Plot Earth and Moon in Earth's rest frame
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:,1,0]*0.0, x[:,1,1]*0.0, x[:,1,2]*0.0)
ax.plot(x[:,2,0]-x[:,1,0], x[:,2,1]-x[:,1,1], x[:,2,2]-x[:,1,2])
plt.show()

