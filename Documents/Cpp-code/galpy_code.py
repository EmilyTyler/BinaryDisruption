from galpy.potential import MWPotential2014, NFWPotential, evaluateDensities
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
import numpy as np

#o = Orbit(vxvv=[1.,0.1,1.1,0.,0.1,0.])
o = Orbit.from_name('sol')

#End time in Gyr
T = 10.0
#Time step size
dt = 0.01
ts = np.arange(0, T+dt, dt, dtype=float)
N_t = np.size(ts)

#Integrate orbit
o.integrate(ts, MWPotential2014)

#Evaluate average dark matter density
#Dark matter potential
nfw = NFWPotential(amp=1.0, a=1.0, normalize=True)
avg_density = 0.0
for i in range(N_t):
    avg_density += evaluateDensities(nfw, o.R(ts[i]), o.z(ts[i]), o.phi(ts[i]), ts[i])*dt/T
print('Average density =', avg_density*bovy_conversion.dens_in_msolpc3(220.0, 8.0))
