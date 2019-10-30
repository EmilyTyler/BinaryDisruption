#To plot test data
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import au, parsec, giga, year
from matplotlib import animation

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

'''
data = np.zeros(10**7, dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number] = float(row[0])
		row_number += 1

N_bins = 1000

d_min = np.min(data)
d_max = np.max(data)
print('min =', d_min)
print('max =', d_max)
dd = (d_max - d_min)/(N_bins)
d_bins = np.array([d_min + i*dd for i in range(N_bins)])
N_d = np.zeros(N_bins)
for i in range(np.size(data)):
	j = int(np.floor((data[i]-d_min)/dd))
	if (data[i] == d_max):
		j = N_bins-1
	N_d[j] += 1
N_d /= np.size(data)
#Move bins into centre for plotting and calculations
d_bins += 0.5*dd
plt.plot(d_bins/au, N_d*au/dd, label='Simulation')
e=0.3
a = 10.0**5.0 * au
plt.plot(d_bins/au, au*1.0/np.pi*d_bins/a*(a**2.0*e**2.0-(d_bins-a)**2.0)**(-1/2), label='Analytic')
#plt.plot(d_bins, 1.0/(2.0*np.pi)*(1.0-e)**(3.0/2.0)*(1.0+e*np.cos(d_bins))**(-2.0))
plt.xlabel(r'$r$, au')
plt.ylabel(r'Probability density, au$^{-1}$')
plt.legend()
plt.show()




#Vector direction tests
data = np.zeros((5*10**6, 3), dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number, 0] = float(row[0])
		data[row_number, 1] = float(row[1])
		data[row_number, 2] = float(row[2])
		row_number += 1

#Plot 3D
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
for i in range(1000):
	ax1.scatter(data[i,0], data[i,1], data[i,2])
plt.show()

#Test theta = arctan(y/x), should be uniform [0,2pi)
theta = np.arctan2(data[:,1], data[:,0])
N_bins = 100
t_min = 0.0
t_max = 2.0*np.pi
dt = (t_max - t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])
N_t = np.zeros(N_bins)
for i in range(np.size(theta)):
	j = int(np.floor((theta[i]-t_min)/dt))
	N_t[j] += 1
plt.plot(t_bins, N_t/dt/np.size(theta))
plt.plot(t_bins, [1.0/N_bins/dt]*N_bins)
plt.show()
#Test z, should be uniform [-1,1]
z = data[:,2]
N_bins = 100
z_min = -1.0
z_max = 1.1
dz = (z_max - z_min)/(N_bins-1)
z_bins = np.array([z_min + i*dz for i in range(N_bins)])
print('z_bins =', z_bins)
N_z = np.zeros(N_bins)
for i in range(np.size(z)):
	j = int(np.floor((z[i]-z_min)/dz))
	N_z[j] += 1
plt.plot(z_bins, N_z/dz/np.size(z))
plt.plot(z_bins, [1.0/N_bins/dz if (z_bins[i]<=1.0) else 0.0 for i in range(N_bins)])
plt.show()
'''

'''
#2D vdv and i distribution
N_data=10**6
data = np.zeros((N_data, 2), dtype=float)

def load(filename):
	with open(filename) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_number = 0
		for row in reader:
			data[row_number] = np.array([float(row[0]), float(row[1])])
			row_number += 1
	return

load('vdv_i_b10e5au.csv')

N_bins = 100
vdv_min = np.min(data[:,0])
vdv_max = np.max(data[:,0])
dvdv = (vdv_max - vdv_min)/(N_bins-1)
vdv_bins = np.array([vdv_min + i*dvdv for i in range(N_bins)])
#Move to centre of bins for plotting
vdv_bins += 0.5*dvdv
N_vdv = np.zeros(N_bins)
for i in range(N_data):
	j = int(np.floor((data[i,0]-vdv_min)/dvdv))
	N_vdv[j] += 1
plt.plot(vdv_bins, N_vdv)
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$, J')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()

N_bins = 100
i_min = 0.0
i_max = np.pi/2.0
di = (i_max - i_min)/(N_bins-1)
i_bins = np.array([i_min + i*di for i in range(N_bins)])
#Move to centre of bins for plotting
i_bins += 0.5*di
N_i = np.zeros(N_bins)
for i in range(N_data):
	j = int(np.floor((data[i,1]-i_min)/di))
	N_i[j] += 1
plt.plot(i_bins, N_i)
plt.xlabel('Perturber inclination, rad')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()

#Contour plot
N_bins = 100
N_vdvi = np.zeros((N_bins, N_bins))
for i in range(np.size(data[:,0])):
	j = int(np.floor((data[i,0]-vdv_min)/dvdv))
	k = int(np.floor((data[i,1]-i_min)/di))
	N_vdvi[j,k] += 1
levels = np.linspace(np.min(N_vdvi), np.max(N_vdvi), num=25)
#N_dv /= np.size(data[:,1])
cs = plt.contourf(vdv_bins, i_bins, np.transpose(N_vdvi), levels=levels)
cbar = plt.colorbar(cs)
cbar.ax.set_ylabel('Number of encounters')
ax = plt.gca()
ax.set_yticks([0., .25*np.pi, .5*np.pi])
ax.set_yticklabels(["$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$"])
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$, J')
plt.ylabel(r'Perturber inclination $i$, rad')
plt.show()


x1 = np.zeros((0,3))
x2 = np.zeros((0,3))

with open('test_nbody.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		x1 = np.append(x1, [[float(row[0]), float(row[1]), float(row[2])]], axis=0)
		x2 = np.append(x2, [[float(row[3]), float(row[4]), float(row[5])]], axis=0)
		row_number += 1
		#if (row_number>100000):
		#	break

#Plot 3D
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
#ax1.plot3D(x1[:,0], x1[:,1], x1[:,2])
#ax1.plot3D(x2[:,0], x2[:,1], x2[:,2])
ax1.plot3D(x1[:,0]-x2[:,0], x1[:,1]-x2[:,1], x1[:,2]-x2[:,2])
plt.show()
'''
'''
N_data = 15018
data = np.zeros(N_data)

with open('dE_break_binary_1000Msol.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number] = float(row[0])
		row_number += 1

N_bins = 10000
d_min = np.min(data)
d_max = np.max(data)
print('min =', d_min)
print('max =', d_max)
dd = (d_max - d_min)/(N_bins)
d_bins = np.array([d_min + i*dd for i in range(N_bins)])
N_d = np.zeros(N_bins)
for i in range(np.size(data)):
	j = int(np.floor((data[i]-d_min)/dd))
	if (data[i] == d_max):
		j = N_bins-1
	N_d[j] += 1
#N_d /= np.size(data)
#Move bins into centre for plotting and calculations
d_bins += 0.5*dd
plt.plot(d_bins, N_d)
plt.xlabel(r'$\Delta E$, au')
plt.ylabel(r'Number')
plt.legend()
plt.show()
'''

'''
N_bins = 100
a_min = 10.0**(-5.0) * 0.17*parsec
a_max = 10.0**4.0 * 0.17*parsec
da = (np.log(a_max) - np.log(a_min))/(N_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_bins)])
N_a = np.zeros(N_bins)

def loadData(file_name, plot_label):
	with open(file_name) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		for row in reader:
			a = float(row[2])
			e = float(row[3])
			if ((a >= 0.0) and (e < 1.0)):
				if (a == a_max):
					j = N_bins-1
				else:
					j = int(np.floor((np.log(a/a_min))/da))
				N_a[j] += 1
	plt.plot((a_bins+0.5*da)/(0.17*parsec), N_a, label = plot_label)

#loadData('ionised_code_a_dist.csv', 'Hyperbolic code')
#loadData('xv_code_a_dist.csv', 'Integration code')
loadData('ionised_code_a_dist_Nbin10e5.csv', 'Hyperbolic code')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel(r'$a/r_J$')
plt.ylabel('Number of binaries')
plt.legend()
plt.show()
'''

'''
#Plot a and e over time
a = np.zeros(0)
r = np.zeros(0)
e = np.zeros(0)
t = np.zeros(0)

with open('binary_rebound_at_990pc.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		a = np.append(a, float(row[0]))
		r = np.append(r, float(row[1]))
		e = np.append(e, float(row[2]))
		t = np.append(t, float(row[3]))


fig, ax1 = plt.subplots()

color = 'red'
ax1.set_xlabel('Time, Gyr')
ax1.set_ylabel('Parsec', color=color)
ax1.plot(t/(giga*year), a/parsec, color=color, label='Semi-major axis')
#ax1.plot(t/(giga*year), r/parsec, color=color, linestyle='--', label = 'Separation')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'dodgerblue'
ax2.set_ylabel('Eccentricity', color=color)
ax2.plot(t/(giga*year), e, color=color)
ax2.tick_params(axis='y', labelcolor=color)

#plt.legend()
fig.tight_layout()
plt.show()
'''

'''
#Plot distribution of unbound binary separations
r_10 = np.zeros(85068, dtype=float)
with open('final_seps_unbound_binaries_10Msol.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		r_10[row_number] = float(row[0])
		row_number += 1

r_1 = np.zeros(3932, dtype=float)
with open('final_seps_unbound_binaries_1Msol.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		r_1[row_number] = float(row[0])
		row_number += 1

N_bins = 10



r_min = 0.4*parsec
r_max = 1000.0*parsec
print('min =', r_min/parsec)
print('max =', r_max/parsec)
dr = (np.log(r_max) - np.log(r_min))/(N_bins)
r_bins = np.array([r_min*np.exp(i*dr) for i in range(N_bins)])
N_r1 = np.zeros(N_bins)
N_r10 = np.zeros(N_bins)
for i in range(np.size(r_1)):
	j = int(np.floor(np.log(r_1[i]/r_min)/dr))
	if (r_1[i] < r_max):
		N_r1[j] += 1
for i in range(np.size(r_10)):
	j = int(np.floor(np.log(r_10[i]/r_min)/dr))
	if (r_10[i] < r_max):
		N_r10[j] += 1

print(N_r10)

N_r1 /= np.sum(N_r1)
N_r10 /= np.sum(N_r10)


#Move bins into centre for plotting and calculations
r_bins += 0.5*dr
plt.semilogx(r_bins/parsec, N_r1, label=r'$M_p=1M_\odot$')
plt.semilogx(r_bins/parsec, N_r10, label=r'$M_p=10M_\odot$')
plt.xlabel(r'Separation, pc')
plt.xlim([r_min/parsec, r_max/parsec])
plt.ylabel(r'Number of binaries')
plt.legend()
plt.show()
'''
'''
#Plot distribution of unbound binary separations over time
t_min = 0.0
t_max = 10.0*giga*year
N_t_bins = 100
dt = (t_max - t_min)/N_t_bins
t_bins = np.array([t_min + i*dt for i in range(N_t_bins)])

r_min = 0.4*parsec
r_max = 1000.0*parsec
N_r_bins = 500
dr = (np.log(r_max) - np.log(r_min))/(N_r_bins)
r_bins = np.array([r_min*np.exp(i*dr) for i in range(N_r_bins)])

N_r1 = np.zeros((N_t_bins, N_r_bins), dtype=float)
N_r10 = np.zeros((N_t_bins, N_r_bins), dtype=float)
N_r100 = np.zeros((N_t_bins, N_r_bins), dtype=float)

t_min_actual=10.0*giga*year
binary_number_previous = -1
i_previous = -1
with open('final_seps_unbound_binaries_1Msol_with_t_10e4bin.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		if (np.size(row)>0):
			#print(row)
			r = float(row[0])
			t = float(row[1])
			binary_number = int(row[2])
			#print(binary_number)
			i = int(np.floor((t - t_min)/dt))
			#print("binary number = ", binary_number)
			#print("i =", i)
			#print("j =", j)
			t_min_actual=np.min([t_min_actual, t])
			if (t == t_max):
				i = N_t_bins-1
			elif (t>t_max):
				continue
			if ((binary_number != binary_number_previous) or (i!=i_previous)):
				#print(row)
				j = int(np.floor(np.log(r/r_min)/dr))
				if (r == r_max):
					j = N_r_bins-1
				N_r1[i,j] += 1
				#print("Increment")
				binary_number_previous = binary_number
				i_previous = i
			#input()
#print(N_r1[N_t_bins-3])
#print(N_r1[N_t_bins-2])
print('t_min_actual =', t_min_actual/(giga*year))


t_min_actual=10.0*giga*year
binary_number_previous = -1
i_previous = -1
with open('final_seps_unbound_binaries_10Msol_with_t_10e4bin.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		if (np.size(row)>0):
			r = float(row[0])
			t = float(row[1])
			binary_number = int(row[2])
			#print(binary_number)
			i = int(np.floor((t - t_min)/dt))	
			t_min_actual=np.min([t_min_actual, t])	
			if (t == t_max):
				i = N_t_bins-1
			elif (t>t_max):
				continue
			if ((binary_number != binary_number_previous) or (i!=i_previous)):
				j = int(np.floor(np.log(r/r_min)/dr))
				if (r == r_max):
					j = N_r_bins-1
				#print('Increment')
				N_r10[i,j] += 1
				binary_number_previous = binary_number
				i_previous = i
				#input()
#print(N_r10[N_t_bins-1])]
print('t_min_actual =', t_min_actual/(giga*year))


t_min_actual=10.0*giga*year
binary_number_previous = -1
i_previous = -1
with open('final_seps_unbound_binaries_100Msol_with_t_10e4bin.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		if (np.size(row)>0):
			r = float(row[0])
			t = float(row[1])
			binary_number = int(row[2])
			i = int(np.floor((t - t_min)/dt))
			t_min_actual=np.min([t_min_actual, t])
			if (t == t_max):
				i = N_t_bins-1
			elif (t>t_max):
				continue
			if ((binary_number != binary_number_previous) or (i!=i_previous)):
				j = int(np.floor(np.log(r/r_min)/dr))
				if (r == r_max):
					j = N_r_bins-1
				N_r100[i,j] += 1
				binary_number_previous = binary_number
				i_previous = i
print('t_min_actual =', t_min_actual/(giga*year))

#Normalise
for i in range(N_t_bins):
	N_r1[i] /= max([1, np.sum(N_r1[i])])
	N_r10[i] /= max([1, np.sum(N_r10[i])])
	N_r100[i] /= max([1, np.sum(N_r100[i])])
#Move bins into centre for plotting and calculations
r_bins += 0.5*dr
t_bins += 0.5*dt

y_max = 0.05

plt.semilogx(r_bins/parsec, N_r1[N_t_bins-1], label=r'$M_p=1M_\odot$')
plt.semilogx(r_bins/parsec, N_r10[N_t_bins-1], label=r'$M_p=10M_\odot$')
plt.semilogx(r_bins/parsec, N_r100[N_t_bins-1], label=r'$M_p=100M_\odot$')
plt.xlabel(r'Separation, pc')
plt.xlim([r_min/parsec, r_max/parsec])
plt.ylim([-0.005,y_max])
plt.ylabel(r'Fraction of currently broken binaries')
plt.legend()
plt.show()

#Generate animation
base_interval = 200
fig = plt.figure()
ax=plt.gca()
plt.xlim([r_min/parsec, r_max/parsec])
plt.ylim([-0.005,y_max])
plt.xlabel(r'Separation, pc')
plt.ylabel(r'Fraction of currently broken binaries')
plt.text(0.05, 0.9, r'$t = {}$Gyr'.format(round(t_bins[0]/(giga*year), 1)), transform=ax.transAxes)
graph = ax.semilogx(r_bins/parsec, N_r1[0], label=r'$M_p=1M_\odot$', color = 'dodgerblue')
graph = ax.semilogx(r_bins/parsec, N_r10[0], label=r'$M_p=10M_\odot$', color = 'darkorange')
graph = ax.semilogx(r_bins/parsec, N_r100[0], label=r'$M_p=100M_\odot$', color = 'forestgreen')
plt.legend(loc='upper right')
def update(i):
	ax.cla()
	plt.xlim([r_min/parsec, r_max/parsec])
	plt.ylim([-0.005,y_max])
	plt.xlabel(r'Separation, pc')
	plt.ylabel(r'Fraction of currently broken binaries')
	plt.text(0.05, 0.9, r'$t = {}$Gyr'.format(round(t_bins[i]/(giga*year), 3)), transform=ax.transAxes)
	graph = ax.semilogx(r_bins/parsec, N_r1[i], label=r'$M_p=1M_\odot$', color = 'dodgerblue')
	graph = ax.semilogx(r_bins/parsec, N_r10[i], label=r'$M_p=10M_\odot$', color = 'darkorange')
	graph = ax.semilogx(r_bins/parsec, N_r100[i], label=r'$M_p=100M_\odot$', color = 'forestgreen')
	plt.legend(loc='upper right')
	return(graph)        
anim = animation.FuncAnimation(fig, update, frames = range(0, N_t_bins), interval=base_interval, repeat=True, repeat_delay=600)
anim.save('unbound_distribution.mp4', writer=writer)
plt.show()
'''


#Plot Jiang and Tremaine fig 4
t_min = 0.0
t_max = 10.0*giga*year
N_t_bins = 100
dt = (t_max - t_min)/N_t_bins
t_bins = np.array([t_min + i*dt for i in range(N_t_bins)])

r_min = 0.4*parsec
r_max = 1000.0*parsec
N_r_bins = 500
dr = (np.log(r_max) - np.log(r_min))/(N_r_bins)
r_bins = np.array([r_min*np.exp(i*dr) for i in range(N_r_bins)])
dr_log = np.zeros(N_r_bins)
for i in range(N_r_bins-1):
	dr_log[i] = r_bins[i+1]-r_bins[i]
dr_log[N_r_bins-1] = r_max - r_bins[N_r_bins-1]

N_r1 = np.zeros((N_t_bins, N_r_bins), dtype=float)

binary_number_previous = -1
i_previous = -1
with open('final_seps_unbound_binaries_1Msol_with_t_10e4bin.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		if (np.size(row)>0):

			r = float(row[0])
			t = float(row[1])
			binary_number = int(row[2])
			i = int(np.floor((t - t_min)/dt))

			if (t == t_max):
				i = N_t_bins-1
			elif (t>t_max):
				continue
			if ((binary_number != binary_number_previous) or (i!=i_previous)):
				j = int(np.floor(np.log(r/r_min)/dr))
				if (r == r_max):
					j = N_r_bins-1
				N_r1[i,j] += 1
				binary_number_previous = binary_number
				i_previous = i




#Normalise
for i in range(N_t_bins):
	N_r1[i] /= max([1, np.sum(N_r1[i])])*(dr_log/(parsec))

#Move bins into centre for plotting and calculations
r_bins += 0.5*dr
t_bins += 0.5*dt

y_max = 0.05

plt.plot(np.log10(r_bins/(1.7*parsec)), N_r1[N_t_bins-1], label=r'$a_i=0.59r_J$')
plt.xlabel(r'$\mathrm{log}_{10}(r/r_J)$')
plt.xlim([-2, 4])
plt.ylim([0,0.8])
plt.ylabel(r'Probability Density')
plt.legend()
plt.show()