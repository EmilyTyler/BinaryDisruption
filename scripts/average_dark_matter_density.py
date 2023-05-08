# Calculates and plots the average dark matter density for binaries in the AMR catalog

import csv
import scipy
from galpy.potential import MWPotential2014, evaluateDensities
from galpy.potential.mwpotentials import DehnenBinney98I
from galpy.orbit import Orbit
from galpy.util import conversion
from galpy.util.conversion import get_physical
from astropy import units
from astropy.units import imperial
imperial.enable()
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=14)

#Integration time
ts = np.linspace(0.,-10.0,3001)*units.Gyr
N_t = np.size(ts)

#Initialise orbits
print('Importing data')
#Import labels
labels = np.zeros(0)
with open("data/AMR_data_names.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        labels = np.append(labels, row[0])
#Find which stars have orbit data
N_orbits = 0
for i in range(np.size(labels)):
    try:
        Orbit.from_name(labels[i], **get_physical(MWPotential2014))
        N_orbits += 1
    except Exception as e:
        print("The error is: ",e)
        print('No orbit for:', labels[i])
        labels[i] = 'no_orbit'
labels = np.where(labels=='NLTT 7795', 'no_orbit', labels)
N_orbits -= 1
labels = np.delete(labels, np.nonzero(labels=='no_orbit'))
#Import data
orbits = [Orbit.from_name(labels[i], **get_physical(MWPotential2014)) for i in range(N_orbits)]
print('Number of binaries with orbits =', N_orbits)

#Dark matter density at sol
dm_density_at_sol = evaluateDensities(MWPotential2014[2], 1.0, 0.0)*conversion.dens_in_msolpc3(**get_physical(MWPotential2014))
print('Dark matter density at solar position =', dm_density_at_sol)

#Integrate orbits
print('Integrating orbits')
for orbi in orbits:
    orbi.turn_physical_off()
    orbi.integrate(ts, MWPotential2014)


plt.figure(figsize=(7.5,5))
#Evaluate time-averaged dark matter density and stellar density
print('Calculating average densities')
avg_dm_density = np.zeros(N_orbits)
avg_stellar_density = np.zeros(N_orbits)
fraction_of_lifetime_in_disk = np.zeros(N_orbits)
Rs = np.zeros(N_orbits)
zs = np.zeros(N_orbits)
for orbi,label,j in zip(orbits,labels, range(N_orbits)):
    #print('Orbit', j+1, 'of', N_orbits)
    orbi.turn_physical_on(**get_physical(MWPotential2014))
    Rs[j] = orbi.rap()/units.kpc
    zs[j] = orbi.zmax()/units.kpc

    for orbi2,label2 in zip(orbits,labels):
        if ((orbi.rap()/units.kpc == orbi2.rap()/units.kpc) and (label!=label2)):
            print('Has duplicate:', label, 'and', label2)

    for i in range(N_t):
        if (abs(orbi.z(ts[i].value)/units.kpc)<0.5):
            fraction_of_lifetime_in_disk[j] += 1/N_t

    for i in range(N_t):
        avg_dm_density[j] += evaluateDensities(MWPotential2014[2], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*conversion.dens_in_msolpc3(220.0, 8.0)
        avg_stellar_density[j] += evaluateDensities(MWPotential2014[0], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*conversion.dens_in_msolpc3(220.0, 8.0) + evaluateDensities(MWPotential2014[1], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*conversion.dens_in_msolpc3(220.0, 8.0)
    
#Add in Quinn et al. initial conditions NLTT 15501
orbit15501 = Orbit([85.91593, 49.38367, 0.210, 81, -176, 262.3], radec=True, **get_physical(DehnenBinney98I))
orbit15501.turn_physical_off()
orbit15501.integrate(ts, MWPotential2014)
orbit15501.turn_physical_on(**get_physical(MWPotential2014))
orbit15501.plot(d1='R', d2='z', **get_physical(MWPotential2014), overplot=True, label='NLTT 15501*', lw=0.4)
plt.xlabel(r'$R$/kpc')
plt.ylabel(r'$z$/kpc')
plt.legend()
plt.tight_layout()
plt.show()

print('Median apocentre distance, kpc =', np.median(Rs))
print('Median maximum z, kpc =', np.median(zs))

print('Highest average dm density:', labels[np.argmax(avg_dm_density)], np.max(avg_dm_density))
print('Lowest average dm density:', labels[np.argmin(avg_dm_density)], np.min(avg_dm_density))

print('Highest average stellar density:', labels[np.argmax(avg_stellar_density)], np.max(avg_stellar_density))
print('Lowest average stellar density:', labels[np.argmin(avg_stellar_density)], np.min(avg_stellar_density))

#Make distribution of average densities
print('Plotting distribution of dark matter densities')
plt.figure(figsize=(7.5,5))
#Setting up bins
rho_min = min(np.min(avg_dm_density), np.min(avg_stellar_density))
rho_max = max(np.max(avg_dm_density), np.max(avg_stellar_density))
N_bins = 270
d_rho = (rho_max-rho_min)/(N_bins-1)
rho_bins = np.array([rho_min + i*d_rho for i in range(N_bins)])
#Bin densities
N_rho_dm = np.zeros(N_bins)
N_rho_stellar = np.zeros(N_bins)
for i in range(N_orbits):
    j = int(np.floor((avg_dm_density[i]-rho_min)/d_rho))
    N_rho_dm[j] += 1
    k = int(np.floor((avg_stellar_density[i]-rho_min)/d_rho))
    N_rho_stellar[k] += 1
N_rho_dm /= N_orbits
N_rho_stellar /= N_orbits
N_rho_dm /= d_rho
N_rho_stellar /= d_rho
#Spline for FWHM
HM=np.max(N_rho_dm)/2
spline = scipy.interpolate.UnivariateSpline(rho_bins, N_rho_dm-HM)
plt.plot(rho_bins, N_rho_dm, label='DM')
plt.plot(rho_bins, N_rho_stellar, label='Stellar')
plt.plot([dm_density_at_sol]*N_bins, N_rho_dm, label='Density at solar radius')
plt.xlabel(r'Time-averaged density / $M_\odot$pc$^{-3}$')
plt.ylabel('Probability density')
plt.legend()
plt.tight_layout()
plt.show()

print('Mean of average dm densities =', np.mean(avg_dm_density))
print('Median of average dm densities =', np.median(avg_dm_density))
print('Mode of average dm densities =', rho_bins[np.argmax(N_rho_dm)], 'to', rho_bins[np.argmax(N_rho_dm)+1])

print('Mean of average densities =', np.mean(avg_stellar_density))
print('Median of average densities =', np.median(avg_stellar_density))
print('Mode of average densities =', rho_bins[np.argmax(N_rho_stellar)], 'to', rho_bins[np.argmax(N_rho_stellar)+1])

print('Maximum apocentre distance, kpc =', np.max(np.sqrt(Rs**2.0 + zs**2.0)))

print('Root mean square of density distribution:', np.sqrt(np.sum(N_rho_dm**2.0)/np.size(N_rho_dm)))
r1, r2 = spline.roots()
print('FWHM of DM density distribution:', r2-r1)

#Plot of avg_dm_density/solar_dm_density vs fraction of time in disk
plt.figure(figsize=(7.5,5))
plt.scatter(fraction_of_lifetime_in_disk, avg_dm_density/dm_density_at_sol, marker='x')
plt.xlabel(r'Fraction of time spent with $|z| < 500$pc')
plt.ylabel(r'Time-averaged dark matter density / $\rho_\odot$')
plt.tight_layout()
plt.show()

print('Finished')