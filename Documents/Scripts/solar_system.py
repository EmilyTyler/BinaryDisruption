import numpy as np
import os
os.system("python setup.py build_ext --inplace")


from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from scipy.constants import year, au, day, G, kilo
from evolve_binary import integrateBinary
import itertools as it
from orbital_elements import orbitalElements

print('Initialising')
#Time array
t = np.array([0.0])
#Finish time
t_end = kilo * year
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
x[:,9:] = x[:,9:]/day

#Masses in kg
m = np.array([1.988544*10.0**30.0, 3.302*10.0**23.0, 48.685*10.0**23.0, 5.97219*10.0**24.0, 6.4185*10.0**23.0, 1898.13*10.0**24.0, 5.68319*10.0**26.0, 86.8103*10.0**24.0, 102.41*10.0**24.0])

print('Integrating orbits')
#Initialise counter
i = 1
while t[i-1] < t_end:
        (x_new, dt) = integrateBinary(9, x[i-1], m, n=20)
        x = np.append(x, [x_new], axis=0)
        t = np.append(t, t[i-1]+dt)
        #Increment counter
        i += 1

#Plot paths in sun's rest frame
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.view_init(elev=90.0)
ax.plot(x[:,0,0]*0.0, x[:,0,1]*0.0, x[:,0,2]*0.0, label='Sun')
ax.plot(x[:,1,0]-x[:,0,0], x[:,1,1]-x[:,0,1], x[:,1,2]-x[:,0,2], label='Mercury')
ax.plot(x[:,2,0]-x[:,0,0], x[:,2,1]-x[:,0,1], x[:,2,2]-x[:,0,2], label='Venus')
ax.plot(x[:,3,0]-x[:,0,0], x[:,3,1]-x[:,0,1], x[:,3,2]-x[:,0,2], label='Earth')
ax.plot(x[:,4,0]-x[:,0,0], x[:,4,1]-x[:,0,1], x[:,4,2]-x[:,0,2], label='Mars')
ax.plot(x[:,5,0]-x[:,0,0], x[:,5,1]-x[:,0,1], x[:,5,2]-x[:,0,2], label='Jupiter')
ax.plot(x[:,6,0]-x[:,0,0], x[:,6,1]-x[:,0,1], x[:,6,2]-x[:,0,2], label='Saturn')
ax.plot(x[:,7,0]-x[:,0,0], x[:,7,1]-x[:,0,1], x[:,7,2]-x[:,0,2], label='Uranus')
ax.plot(x[:,8,0]-x[:,0,0], x[:,8,1]-x[:,0,1], x[:,8,2]-x[:,0,2], label='Neptune')
ax.set_xlim3d([-50.0*au, 50.0*au])
ax.set_ylim3d([-50.0*au, 50.0*au])
ax.set_zlim([-10.0*au, 10.0*au])
plt.legend()
plt.show()



#Generate 3D animation
base_interval = 1
fig = plt.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=90.0)
ax.set_xlim3d([-30.0*au, 30.0*au])
ax.set_ylim3d([-30.0*au, 30.0*au])
ax.set_zlim([-1.0*au, 1.0*au])
colours = ["yellow", "darkgrey", "gold", "royalblue", "r", "sienna", "orange", "deepskyblue", "slateblue"]
graph = ax.scatter(x[0,:9,0]-x[0,0,0], x[0,:9,1]-x[0,0,1], x[0,:9,2]-x[0,0,2], c=colours)
def update(i):
        ax.clear()
        ax.scatter(x[i,:9,0]-x[i,0,0], x[i,:9,1]-x[i,0,1], x[i,:9,2]-x[i,0,2], c=colours)
        ax.view_init(elev=90.0)
        ax.set_xlim3d([-30.0*au, 30.0*au])
        ax.set_ylim3d([-30.0*au, 30.0*au])
        ax.set_zlim([-1.0*au, 1.0*au])
        #graph = ax.scatter(x[i,:9,0], x[i,:9,1], x[i,:9,2], c=colours)
        #if i < np.size(x[:,0,0]):
                #anim.event_source.interval = base_interval * (t[i+1]-t[i])/t[1]
        return(graph)        
anim = animation.FuncAnimation(fig, update, frames=np.arange(np.size(x[:,0,0])), interval=base_interval, repeat=True)
plt.show()
#anim.save('solar_system.mp4')


#TESTS
N_t = np.size(x[:,0,0])
#Energy conservation
print('Calculating total energy')
E = np.zeros(N_t)
for i in range(N_t):
        for j in range(9):
                #Kinetic energy
                E[i] += 0.5*m[j]*np.dot(x[i,j+9], x[i,j+9])
                #Gravitational energy
                for k in it.chain(range(j), range(j+1, 9)):
                        E[i] += 0.5 * - G*m[j]*m[k]/np.linalg.norm(x[i,j]-x[i,k])
plt.title('Solar System energy conservation test')
plt.xlabel('Time, yr')
plt.ylabel('Fractional energy change')
plt.plot(t/year, (E-E[0])/E[0])
plt.show()

#Angular momentum conservation
print('Calculating total angular momentum')
L = np.zeros((N_t,3))
for i in range(N_t):
        for j in range(9):
                L[i] += m[j]*np.cross(x[i,j], x[i,j+9])
plt.title(r'Solar system angular momentum conservation test: $x$-component')
plt.xlabel('Time, yr')
plt.ylabel('Fractional angular momentum change')
plt.plot(t/year, (L[:,0]-L[0,0])/L[0,0])
plt.show()
plt.title(r'Solar system angular momentum conservation test: $y$-component')
plt.xlabel('Time, yr')
plt.ylabel('Fractional angular momentum change')
plt.plot(t/year, (L[:,1]-L[0,1])/L[0,1])
plt.show()
plt.title(r'Solar system angular momentum conservation test: $z$-component')
plt.xlabel('Time, yr')
plt.ylabel('Fractional angular momentum change')
plt.plot(t/year, (L[:,2]-L[0,2])/L[0,2])
plt.show()

#Semi-major axes and eccentricities of Jupiter and Saturn
print('Calculating a and e for Jupiter and Saturn')
a_Jup = np.zeros(N_t)
a_Sat = np.zeros(N_t)
e_Jup = np.zeros(N_t)
e_Sat = np.zeros(N_t)
for i in range(N_t):
        notBound_Jup, a_Jup[i], e_Jup[i] = orbitalElements(np.array([x[i,5], x[i,0], x[i,14], x[i,9]]), m[5], m[0])
        notBound_Sat, a_Sat[i], e_Sat[i] = orbitalElements(np.array([x[i,6], x[i,0], x[i,15], x[i,9]]), m[6], m[0])
#Plot semi-major axes
fig = plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
plt.title('Semi-major axis test for Jupiter and Saturn')
ax1.set_xlabel('Time, yr')
ax1.set_ylabel('Semi-major axis of Jupiter, au')
ax1.plot(t/year, a_Jup/au, label='Jupiter', color='dodgerblue')
ax2.set_ylabel('Semi-major axis of Saturn, au')
ax2.plot(t/year, a_Sat/au, label='Saturn', color='darkorange')
plt.legend()
plt.show()
#Plot eccentricity
fig = plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
plt.title('Eccentricity test for Jupiter and Saturn')
ax1.set_xlabel('Time, yr')
ax1.set_ylabel('Eccentricity of Jupiter')
ax1.plot(t/year, e_Jup, label='Jupiter', color='dodgerblue')
ax2.set_ylabel('Eccentricity of Saturn')
ax2.plot(t/year, e_Sat, label='Saturn', color='darkorange')
plt.legend
plt.show()

#Eccentricity of Earth to test for Milankovitch cycles
print('Calculating eccentricity of Earth')
e_Earth = np.zeros(N_t)
for i in range(N_t):
        notBound_Earth, a_Earth, e_Earth[i] = orbitalElements(np.array([x[i,3], x[i,0], x[i,12], x[i,9]]), m[3], m[0])
plt.title('Milankovitch test')
plt.xlabel('Time, yr')
plt.ylabel('Eccentricity')
plt.plot(t/year, e_Earth)
plt.show()








print('Finished')




