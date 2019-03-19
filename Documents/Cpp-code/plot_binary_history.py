import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import parsec, giga, year

N_bin = 10
N_enc_max = 10000000

'''
#Save data from csvs
t = np.zeros((N_bin, N_enc_max))
v_rel = np.zeros((N_bin, N_enc_max))
b = np.zeros((N_bin, N_enc_max))
a = np.zeros((N_bin, N_enc_max))
e = np.zeros((N_bin, N_enc_max))
N_enc = np.zeros((N_bin, N_enc_max))
N_enc_max_actual = 0
for i in range(N_bin):
	with open('binary_history_{}.csv'.format(i)) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_number = 0
		for row in reader:
			t[i, row_number] = float(row[0])
			v_rel[i, row_number] = float(row[1])
			b[i, row_number] = float(row[2])
			a[i, row_number] = float(row[3])
			e[i, row_number] = float(row[4])
			N_enc[i, row_number] = float(row[5])
			N_enc_max_actual = np.max([N_enc_max_actual, float(row[5])])
			row_number += 1
print('Max number of encounters =', N_enc_max_actual)
np.savez('binary_history.npz', t=t, v_rel=v_rel, b=b, a=a, e=e, N_enc=N_enc)
'''


#Load data
t = np.load('binary_history.npz')['t']
v_rel = np.load('binary_history.npz')['v_rel']
b = np.load('binary_history.npz')['b']
a = np.load('binary_history.npz')['a']
e = np.load('binary_history.npz')['e']
N_enc = np.load('binary_history.npz')['N_enc']


for i in range(N_bin):
	plt.plot(t[i][np.where(a[i]>0.0)]/(giga*year), a[i][np.where(a[i]>0.0)]/parsec)
	plt.xlabel("Time, Gyr")
	plt.ylabel("Semi-major axis, pc")
	plt.show()

