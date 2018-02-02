import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import year, au, day
from evolve_binary import integrateBinary

#Time array
t = np.array([0.0])
#Finish time
t_end = 1.0 * year
#Initial positions and velocities
#Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
#From http://ssd.jpl.nasa.gov/horizons.cgi#top
#Positions in au
#Velocities in au/day
x = np.array([[[1.611753050316713*10.0**(-3.0), 6.287395902904349*10.0**(-3.0), -1.155915184874496*10.0**(-4.0)],
              [6.714439290556447*10.0**(-2.0), -4.458017597020139*10.0**(-1.0), -4.306906590105406*10.0**(-2.0)],
              [6.058577775513969*10.0**(-1.0), -3.978890766518131*10.0**(-1.0), -4.052872656395953*10.0**(-2.0)],
              [-6.683627241220385*10.0**(-1.0), 7.289407814811039*10.0**(-1.0), -1.442082957601778*10.0**(-4.0)],
              [-1.404071460895274, -7.609046142669975*10.0**(-1.0), 1.830610097198241*10.0**(-2.0)],
              [-4.108739518612015, -3.536036469456197, 1.065687757183578*10.0**(-1.0)],
              [2.165954861587580*10.0**(-1.0), -1.005524322932917*10.0, 1.662086535110869*10.0**(-1.0)],
              [1.766472239120164*10.0, 9.169019219652498, -1.947949360630087*10.0**(-1.0)],
              [2.870962163943246*10.0, -8.494786031095629, -4.867085258202105*10.0**(-1.0)],
              [-6.080295207354824*10.0**(-6.0), 4.786160746927203*10.0**(-6.0), 1.465434120388129*10.0**(-7.0)],
              [2.219454879577484*10.0**(-2.0), 5.484607005960727*10.0**(-3.0), -1.588758482843254*10.0**(-3.0)],
              [1.111478303587297*10.0**(-2.0), 1.673538595552732*10.0**(-2.0), -4.120991197370340*10.0**(-4.0)],
              [-1.289845466177702*10.0**(-2.0), -1.174976924421575*10.0**(-2.0), 1.826214029279664*10.0**(-8.0)],
              [7.221580058529674*10.0**(-3.0), -1.108282156848348*10.0**(-2.0), -4.095639230902802*10.0**(-4.0)],
              [4.833921698442846*10.0**(-3.0), -5.361343981386876*10.0**(-3.0), -8.587955336416548*10.0**(-5.0)],
              [5.270307743004450*10.0**(-3.0), 1.024082636680083*10.0**(-4.0), -2.116232593192775*10.0**(-4.0)],
              [-1.840590014082311*10.0**(-3.0), 3.307483512995010*10.0**(-3.0), 3.618814075040750*10.0**(-5.0)],
              [8.693395514943853*10.0**(-4.0), 3.028723832678445*10.0**(-3.0), -8.212989799798345*10.0**(-5.0)]]])
#Convert to si
x = x*au
x[9:] = x[9:]/day

#Masses in kg
#m = np.array([1.988544*10.0**30.0, 3.302*10.0**23.0, 48.685*10.0**23.0, 5.97219*10.0**24.0, 6.4185*10.0**23.0, 1898.13*10.0**24.0, 5.68319*10.0**26.0, 86.8103*10.0**24.0, 102.41*10.0**24.0])
m = np.array([1.988544*10.0**30.0, 3.302*10.0**23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#Initialise counter
i = 1
while t[i-1] < t_end:
        (x_new, dt) = integrateBinary(9, x[i-1], m, n=1)
        x = np.append(x, [x_new], axis=0)
        t = np.append(t, t[i-1]+dt)
        #Increment counter
        i += 1
        
#Plot paths in sun's rest frame
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:,0,0]*0.0, x[:,0,1]*0.0, x[:,0,2]*0.0, label='Sun')
ax.plot(x[:,1,0]-x[:,0,0], x[:,1,1]-x[:,0,1], x[:,1,2]-x[:,0,2], label='Mercury')
#ax.plot(x[:,2,0]-x[:,0,0], x[:,2,1]-x[:,0,1], x[:,2,2]-x[:,0,2], label='Venus')
#ax.plot(x[:,3,0]-x[:,0,0], x[:,3,1]-x[:,0,1], x[:,3,2]-x[:,0,2], label='Earth')
#ax.plot(x[:,4,0]-x[:,0,0], x[:,4,1]-x[:,0,1], x[:,4,2]-x[:,0,2], label='Mars')
#ax.plot(x[:,5,0]-x[:,0,0], x[:,5,1]-x[:,0,1], x[:,5,2]-x[:,0,2], label='Jupiter')
#ax.plot(x[:,6,0]-x[:,0,0], x[:,6,1]-x[:,0,1], x[:,6,2]-x[:,0,2], label='Saturn')
#ax.plot(x[:,7,0]-x[:,0,0], x[:,7,1]-x[:,0,1], x[:,7,2]-x[:,0,2], label='Uranus')
#ax.plot(x[:,8,0]-x[:,0,0], x[:,8,1]-x[:,0,1], x[:,8,2]-x[:,0,2], label='Neptune')
plt.legend()
plt.show()