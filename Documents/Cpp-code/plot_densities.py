import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import parsec
plt.rc('font', family='serif')

def densityNFW(r):
	rho_0 = 0.055
	r_0 = 8.1*1000.0
	return rho_0 * (r/r_0*(1 + r/r_0)**2.0)**(-1.0)

def densityStellar(r, z):
	rho_0 = 0.00081
	q = 0.6
	r_break = 19.0*1000.0
	if (r <= r_break):
		n = 2.4
	else:
		n = -4.8
	return rho_0 * ((r/1000.0)**2.0 - (z/1000.0)**2.0 + (z/1000.0)**2.0/q)**(n/2.0)

r = np.linspace(20.0*1000.0, 60.0*1000.0, num=10)
z = np.linspace(20.0*1000.0, 60.0*1000.0, num=10)


for i in range(10):
	for j in range(10):
		print()
		print(r[i]/1000.0)
		print(z[j]/1000.0)
		print(densityNFW(r[i]))
		print(densityStellar(r[i], z[j]))
		plt.scatter(r[i]/1000.0, densityNFW(r[i]), color='dodgerblue', marker='x')
		plt.scatter(r[i]/1000.0, densityStellar(r[i], z[j]), color='darkorange', marker='x')

plt.xlabel('Galactic Radius, kpc')
plt.ylabel(r'Mass density, $M_\odot$pc$^{-3}$')
plt.show()