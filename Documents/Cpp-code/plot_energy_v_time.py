import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec, giga, year, G
import matplotlib.ticker as mticker


'''
Es = np.zeros((0,0), dtype=float)
ts = np.zeros((0,0), dtype=float)
i_previous = -1
binary_number = -1
row_number = 0
with open('energy_v_time_bound_large_r_100Msol.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		#print(row_number, binary_number)
		E = float(row[0])
		t = float(row[1])
		i = int(row[2])
		if (i != i_previous):
			#print(Es.shape)
			#print(np.array([np.zeros(np.array(Es.shape)[1])]).shape)
			Es = np.append(Es, [np.zeros(np.array(Es.shape)[1])], axis=0)
			ts = np.append(ts, [np.zeros(np.array(ts.shape)[1])], axis=0)
			i_previous = i
			binary_number += 1
			row_number = 0
			#print(Es)
		while (np.array(Es.shape)[1]-1 < row_number):
			Es = np.append(Es, np.transpose([np.zeros(np.array(Es.shape)[0])]), axis=1)
			ts = np.append(ts, np.transpose([np.zeros(np.array(ts.shape)[0])]), axis=1)
			#print(Es.shape)
		Es[binary_number, row_number] = E
		ts[binary_number, row_number] = t
		row_number += 1
		#input()
print('Start plotting?')
input()
m1 = 2.0*10.0**30.0
m2 = m1
for i in range(np.size(Es[:,0])):
	print('Plot', i+1, 'of', np.size(Es[:,0]))
	plt.plot(ts[i][np.where(ts[i]>0.0)]/(giga*year), Es[i][np.where(ts[i]>0.0)])
	plt.plot(ts[i][np.where(ts[i]>0.0)]/(giga*year), [0.0]*np.size(ts[i][np.where(ts[i]>0.0)]), linestyle='--', color='grey')
	plt.xlabel('Time, Gyr')
	plt.ylabel('Energy, J')
	ax1 = plt.gca()
	ax1.yaxis.set_major_locator(plt.MaxNLocator(15))
	ax2 = ax1.twinx()
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_ylabel('Semi-major axis, pc')
	# apply a function formatter
	formatter = mticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(-G*m1*m2/(2.0*x)/parsec))
	ax2.yaxis.set_major_formatter(formatter)
	ax2.yaxis.set_major_locator(plt.MaxNLocator(15))
	plt.show()
	
'''
Es = np.zeros((0,0), dtype=float)
Es_nbody = np.zeros((0,0), dtype=float)
ts = np.zeros((0,0), dtype=float)
i_previous = -1
binary_number = -1
row_number = 0
with open('energy_v_time_nbody_1Msol_eta0_00000002_n10_1bin.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		#print(row_number, binary_number)
		E = float(row[0])
		E_nbody = float(row[1])
		t = float(row[2])
		i = int(row[3])
		if (i != i_previous):
			#print(Es.shape)
			#print(np.array([np.zeros(np.array(Es.shape)[1])]).shape)
			Es = np.append(Es, [np.zeros(np.array(Es.shape)[1])], axis=0)
			Es_nbody = np.append(Es_nbody, [np.zeros(np.array(Es_nbody.shape)[1])], axis=0)
			ts = np.append(ts, [np.zeros(np.array(ts.shape)[1])], axis=0)
			i_previous = i
			binary_number += 1
			row_number = 0
			#print(Es)
		while (np.array(Es.shape)[1]-1 < row_number):
			Es = np.append(Es, np.transpose([np.zeros(np.array(Es.shape)[0])]), axis=1)
			Es_nbody = np.append(Es_nbody, np.transpose([np.zeros(np.array(Es_nbody.shape)[0])]), axis=1)
			ts = np.append(ts, np.transpose([np.zeros(np.array(ts.shape)[0])]), axis=1)
			#print(Es.shape)
		Es[binary_number, row_number] = E
		Es_nbody[binary_number, row_number] = E_nbody
		ts[binary_number, row_number] = t
		row_number += 1
		#input()
print('Start plotting?')
input()
m1 = 2.0*10.0**30.0
m2 = m1
for i in range(np.size(Es[:,0])):
	print('Plot', i+1, 'of', np.size(Es[:,0]))
	plt.plot(ts[i][np.where(ts[i]>0.0)]/(giga*year), Es[i][np.where(ts[i]>0.0)], label='Hyperbolic equations')
	plt.plot(ts[i][np.where(ts[i]>0.0)]/(giga*year), Es_nbody[i][np.where(ts[i]>0.0)],label='2 body')
	plt.plot(ts[i][np.where(ts[i]>0.0)]/(giga*year), [0.0]*np.size(ts[i][np.where(ts[i]>0.0)]), linestyle='--', color='grey')
	plt.xlabel('Time, Gyr')
	plt.ylabel('Energy, J')
	ax1 = plt.gca()
	ax1.yaxis.set_major_locator(plt.MaxNLocator(15))
	ax2 = ax1.twinx()
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_ylabel('Semi-major axis, pc')
	# apply a function formatter
	formatter = mticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(-G*m1*m2/(2.0*x)/parsec))
	ax2.yaxis.set_major_formatter(formatter)
	ax2.yaxis.set_major_locator(plt.MaxNLocator(15))
	plt.legend()
	plt.show()

