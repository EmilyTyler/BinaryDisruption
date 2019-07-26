import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au
import matplotlib.lines as mlines
plt.rc('font', family='serif')

a_frac_avg = np.zeros((4, 10))
bs = np.zeros((4, 10))

def loadData(index, filename):
	with open(filename) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_number = 0
		for row in reader:
			a_frac_avg[index, row_number] = float(row[0])
			bs[index, row_number] = float(row[1])
			row_number += 1

loadData(0, 'impulse_test_a10e5au_Mp1Msol.csv')
loadData(1, 'impulse_test_a10e5au_Mp10Msol.csv')
loadData(2, 'impulse_test_a10e5au_Mp100Msol.csv')
loadData(3, 'impulse_test_a10e5au_Mp1000Msol.csv')



for i in range(4):
	plt.plot(bs[i]/au, abs(a_frac_avg[i]), linestyle='--')
#ax = plt.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.xlabel('Impact parameter,au')
#plt.ylabel('Absolute average fractional error in semi-major axis change')
#plt.legend()
#plt.show()

loadData(0, 'impulse_test_a10e4au_Mp1Msol.csv')
loadData(1, 'impulse_test_a10e4au_Mp10Msol.csv')
loadData(2, 'impulse_test_a10e4au_Mp100Msol.csv')
loadData(3, 'impulse_test_a10e4au_Mp1000Msol.csv')

plt.gca().set_prop_cycle(None)

for i in range(4):
	plt.plot(bs[i]/au, abs(a_frac_avg[i]), label=r'$M_p=10^{}M_\odot$'.format(i))
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Impact parameter,au')
plt.ylabel('Absolute average fractional error in semi-major axis change')

plt.legend()
plt.show()