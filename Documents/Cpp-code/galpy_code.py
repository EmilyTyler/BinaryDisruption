from galpy.potential.mwpotentials import DehnenBinney98I
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
from galpy.util.bovy_conversion import get_physical
from astropy import units
from astropy.units import imperial
imperial.enable()
import numpy as np
import matplotlib.pyplot as plt

#Integration time
ts = np.linspace(0.,-10.0,2001)*units.Gyr

#Initialise orbits
orbits = [Orbit([49.62049, -7.14044, 0.219, 171, -353, 121.6], radec=True, **get_physical(DehnenBinney98I)), Orbit.from_name('NLTT 10536', **get_physical(DehnenBinney98I)), Orbit([85.91593, 49.38367, 0.210, 81, -176, 262.3], radec=True, **get_physical(DehnenBinney98I)), Orbit([94.91613, -30.70087, 0.348, 328, -172, 268.2], radec=True, **get_physical(DehnenBinney98I)), Orbit.from_name('NLTT 16394', **get_physical(DehnenBinney98I)), Orbit.from_name('NLTT 39456', **get_physical(DehnenBinney98I))]

#Integrate orbits
for orbi in orbits:
    orbi.turn_physical_off()
    orbi.integrate(ts, DehnenBinney98I)

#Labels for plot
labels = ['NLTT 10536 Quinn', 'NLTT 10536 SIMBAD', 'NLTT 15501 Quinn', 'NLTT 16394 Quinn', 'NLTT 16394 SIMBAD', 'NLTT 39456 SIMBAD']
#Plot over a previous plot?
overp = [False, True, True, True, True, True]
#Variable to plot on the x-axis
xplots = ['R', 'R', '-R', 'R', 'R', 'R']

#Plot orbits
for orbi,label,over,xplot in zip(orbits,labels,overp,xplots):
    orbi.plot(d1=xplot, d2='z', **get_physical(DehnenBinney98I), overplot=over, label=label, lw=0.4)
#Add legend to plot
plt.legend()
#Show plot
plt.show()