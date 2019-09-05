from galpy.potential import MWPotential2014, evaluateDensities
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
import numpy as np

#o = Orbit(vxvv=[1.,0.1,1.1,0.,0.1,0.])
o = Orbit()
o.turn_physical_off()

#End time in Gyr
T = 7.0 / bovy_conversion.time_in_Gyr(220.0, 8.0)
#Time step size
dt = 0.001 / bovy_conversion.time_in_Gyr(220.0, 8.0)
ts = np.arange(0, T+dt, dt, dtype=float)
N_t = np.size(ts)

#Integrate orbit
o.integrate(ts, MWPotential2014)

#Evaluate average dark matter density and stellar density
avg_dm_density = 0.0
avg_stellar_density = 0.0
for i in range(N_t):
    avg_dm_density += evaluateDensities(MWPotential2014[2], o.R(ts[i]), o.z(ts[i]), o.phi(ts[i]), ts[i])*dt/T
    avg_stellar_density += evaluateDensities(MWPotential2014[0], o.R(ts[i]), o.z(ts[i]), o.phi(ts[i]), ts[i])*dt/T + evaluateDensities(MWPotential2014[1], o.R(ts[i]), o.z(ts[i]), o.phi(ts[i]), ts[i])*dt/T
print('Average dark matter density =', avg_dm_density*bovy_conversion.dens_in_msolpc3(220.0, 8.0))
print('Average stellar density =', avg_stellar_density*bovy_conversion.dens_in_msolpc3(220.0, 8.0))