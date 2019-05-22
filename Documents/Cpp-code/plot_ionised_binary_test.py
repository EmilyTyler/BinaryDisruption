import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import csv
import numpy as np
from scipy.constants import au
plt.rc('font', family='serif')

N_bin = 10**5
M_p = np.array([10.0**0.0, 10.0**1.0, 10.0**2.0, 10.0**3.0, 10.0**4.0])*2.0*10.0**30.0
a = np.array([10.0**3.0, 10.0**4.0, 10.0**5.0, 10.0**6.0]) * au

#Number of binaries bound that were once unbound
number_rebound = np.array([	[0, 0, 116, 1957],
							[1, 210, 19372, 7062],
							[11, 13102, 42055, 15972],
							[22, 8642, 34301, 26832],
							[2, 291, 4571, 50000]])

#Number of binaries unbound and within 100pc
number_close = np.array([	[33, 342, 965, 1560],
							[278, 3128, 2079, 2596],
							[1320, 12448, 5914, 3738],
							[1484, 40753, 37834, 15606],
							[1530, 54019, 89270, 50000]])

#Error due to e=1
error = np.array([	[0, 0, 9, 2],
					[0, 4, 23, 5],
					[0, 3, 8, 1],
					[0, 0, 3, 1],
					[0, 0, 0, 50000]])


#Plot number rebound
markerstyle = ['^', 'x', 'v', '*', '.']
markercolor = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
for i in range(np.size(M_p)):
	plt.scatter(a/au, number_rebound[i], marker=markerstyle[i], color = markercolor[i])
	plt.errorbar(a/au, number_rebound[i], yerr=[error[i], error[i]], linestyle="None", color = markercolor[i])

#Legend entries
legend_entries = [	Line2D([0], [0], marker=markerstyle[0], color='w', markerfacecolor=markercolor[0], label=r'$M_p = 10^0M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[1], color='w', markerfacecolor=markercolor[1], label=r'$M_p = 10^1M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[2], color='w', markerfacecolor=markercolor[2], label=r'$M_p = 10^2M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[3], color='w', markerfacecolor=markercolor[3], label=r'$M_p = 10^3M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[4], color='w', markerfacecolor=markercolor[4], label=r'$M_p = 10^4M_\odot$', markersize=9)]

plt.legend(handles = legend_entries)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Initial semi-major axis, au')
plt.ylabel('Number of bound binaries that were once unbound')
plt.show()


#Plot number within 100pc
markerstyle = ['^', 'x', 'v', '*', '.']
markercolor = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
for i in range(np.size(M_p)):
	plt.scatter(a/au, number_close[i], marker=markerstyle[i], color = markercolor[i])
	plt.errorbar(a/au, number_close[i], yerr=[error[i], error[i]], linestyle="None", color = markercolor[i])

#Legend entries
legend_entries = [	Line2D([0], [0], marker=markerstyle[0], color='w', markerfacecolor=markercolor[0], label=r'$M_p = 10^0M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[1], color='w', markerfacecolor=markercolor[1], label=r'$M_p = 10^1M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[2], color='w', markerfacecolor=markercolor[2], label=r'$M_p = 10^2M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[3], color='w', markerfacecolor=markercolor[3], label=r'$M_p = 10^3M_\odot$', markersize=9),
					Line2D([0], [0], marker=markerstyle[4], color='w', markerfacecolor=markercolor[4], label=r'$M_p = 10^4M_\odot$', markersize=9)]

plt.legend(handles = legend_entries)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Initial semi-major axis, au')
plt.ylabel(r'Number of unbound binaries with a separation $<100$pc')
plt.show()