import csv
from galpy.potential import MWPotential2014, evaluateDensities
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
from galpy.util.bovy_conversion import get_physical
from astropy import units
from astropy.units import imperial
imperial.enable()
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')

#Integration time
ts = np.linspace(0.,-10.0,3001)*units.Gyr
N_t = np.size(ts)

#Initialise orbits
print('Importing data')
#Import labels
labels = np.zeros(0)
with open("AMR_data_names.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        labels = np.append(labels, row[0])
#Find which stars have orbit data
N_orbits = 0
for i in range(np.size(labels)):
    try:
        Orbit.from_name(labels[i], **get_physical(MWPotential2014))
        N_orbits += 1
    except:
        labels[i] = 'no_orbit'
labels = np.where(labels=='NLTT 7795', 'no_orbit', labels)
N_orbits -= 1
labels = np.delete(labels, np.nonzero(labels=='no_orbit'))
#Import data
orbits = [Orbit.from_name(labels[i], **get_physical(MWPotential2014)) for i in range(N_orbits)]
print('Number of binaries with orbits =', N_orbits)

#Dark matter density at sol
dm_density_at_sol = evaluateDensities(MWPotential2014[2], 1.0, 0.0)*bovy_conversion.dens_in_msolpc3(**get_physical(MWPotential2014))
print('Dark matter density at solar position =', dm_density_at_sol)

#Integrate orbits
print('Integrating orbits')
for orbi in orbits:
    orbi.turn_physical_off()
    orbi.integrate(ts, MWPotential2014)


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
        avg_dm_density[j] += evaluateDensities(MWPotential2014[2], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*bovy_conversion.dens_in_msolpc3(220.0, 8.0)
        avg_stellar_density[j] += evaluateDensities(MWPotential2014[0], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*bovy_conversion.dens_in_msolpc3(220.0, 8.0) + evaluateDensities(MWPotential2014[1], orbi.R(ts[i].value), orbi.z(ts[i].value), phi=orbi.phi(ts[i].value), t=ts[i])/N_t*bovy_conversion.dens_in_msolpc3(220.0, 8.0)
    
    #if ((label == 'NLTT 1715') or (label == 'NLTT 10536') or (label == 'NLTT 15501') or (label == 'NLTT 16394') or (label == 'NLTT 39456')):
    #if ((8.1<Rs[j]< 9.9) and (0.576< zs[j]<0.704)):
    #if (zs[j]> 0.5):
    #if (fraction_of_lifetime_in_disk[j] < 0.159):
        #print(label)
        #print('Average dark matter density =', avg_dm_density[j])
        #print('Average stellar density =', avg_stellar_density[j])
        #Plot orbit to make sure it makes sense to calculate average density
        #orbi.plot(d1='R', d2='z', **get_physical(MWPotential2014), overplot=False, label=label, lw=0.4)
        #plt.legend()
        #plt.show()
        #orbi.turn_physical_on(**get_physical(MWPotential2014))
        #orbi.turn_physical_on(ro=8., vo=220.)
        #print('Orbital parameters:')
        #print('Energy, 100km^2 s^-2 =', orbi.E()/100.0)
        #print('Total angular Momentum, 10km kpc s^-1 =', np.sqrt(orbi.L()[0]**2.0 + orbi.L()[1]**2.0 + orbi.L()[2]**2.0)/10.0)
        #print('Pericenter distance, kpc =', orbi.rperi())
        #print('Maximum z, kpc =', orbi.zmax())
        #print('Apocentre distance, kpc =', orbi.rap())
        #print('Fraction of lifetime in disk =', fraction_of_lifetime_in_disk[j])
    
print('Median apocentre distance, kpc =', np.median(Rs))
print('Median maximum z, kpc =', np.median(zs))

#Save average densities data
np.savez('average_densities.npz', avg_dm_density=avg_dm_density, avg_stellar_density=avg_stellar_density)
#Load average densities data
#npzfile = np.load('average_densities.npz')
#avg_dm_density = npzfile['avg_dm_density']
#avg_stellar_density = npzfile['avg_stellar_density']

print('Highest average dm density:', labels[np.argmax(avg_dm_density)], np.max(avg_dm_density))
print('Lowest average dm density:', labels[np.argmin(avg_dm_density)], np.min(avg_dm_density))

print('Highest average stellar density:', labels[np.argmax(avg_stellar_density)], np.max(avg_stellar_density))
print('Lowest average stellar density:', labels[np.argmin(avg_stellar_density)], np.min(avg_stellar_density))

#Make distribution of average densities
print('Plotting distribution of dark matter densities')
#Setting up bins
rho_min = np.min(avg_dm_density)
rho_max = np.max(avg_dm_density)
N_bins = 20
d_rho = (rho_max-rho_min)/(N_bins-1)
rho_bins = np.array([rho_min + i*d_rho for i in range(N_bins)])
#Bin densities
N_rho = np.zeros(N_bins)
for i in range(N_orbits):
    j = int(np.floor((avg_dm_density[i]-rho_min)/d_rho))
    N_rho[j] += 1
N_rho /= N_orbits
N_rho /= d_rho
plt.plot(rho_bins, N_rho)
plt.plot([dm_density_at_sol]*N_bins, N_rho, label='Density at solar radius')
plt.xlabel(r'Time-averaged dark matter density, $M_\odot$pc$^{-3}$')
plt.ylabel('Probability density')
plt.legend()
plt.show()

print('Mean of average dm densities =', np.mean(avg_dm_density))
print('Median of average dm densities =', np.median(avg_dm_density))
print('Mode of average dm densities =', rho_bins[np.argmax(N_rho)], 'to', rho_bins[np.argmax(N_rho)+1])

print('Mean of average densities =', np.mean(avg_stellar_density))
print('Median of average densities =', np.median(avg_stellar_density))
#print('Mode of average densities =', rho_bins[np.argmax(N_rho)], 'to', rho_bins[np.argmax(N_rho)+1])

print('Fraction of lifetimes in disk =', fraction_of_lifetime_in_disk)

'''
plt.scatter(fraction_of_lifetime_in_disk, avg_dm_density/dm_density_at_sol, marker='x')
plt.xlabel(r'Fraction of time spent with $|z| < 500$pc')
plt.ylabel(r'Time-averaged dark matter density, $\rho_\odot$')
plt.show()
'''
'''
#Make file of labels, t_d/t
with open('AMR_data_names_galpyt_d.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for label, fraction in zip(labels, fraction_of_lifetime_in_disk):
        writer.writerow([label, fraction])
'''
print('Finished')