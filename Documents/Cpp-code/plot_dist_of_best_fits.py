import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rc('font', family='serif')

Mps_rproj_dist_random = np.array([[89,89,70,70,50,79,79,89,70,100],
                                  [79,79,79,79,89,100,79,79,79,79],
                                  [100,100,100,100,100,100,79,89,100,100]])
Mps_s_dist_100pc = np.array([[100,70,70,70,70,70,70,70,70,70],
                             [50,50,70,70,70,70,70,70,70,70],
                             [79,89,100,79,100,100,79,79,100,79]])
Mps_s_dist_random = np.array([[1000,300,300,500,1000,840,79,500,840,500],
                              [100,500,707,300,89,500,840,707,707,840],
                              [840,595,707,707,840,1000,707,840,840,500]])

data = np.append(Mps_rproj_dist_random, Mps_s_dist_random, axis=0)
data = np.append(data, Mps_s_dist_100pc, axis=0)
labels = ['Projected separation, 251 binaries', 'Projected separation, 173 binaries', 'Projected separation, 138 binaries',
          'Angular separation, random distance, 251 binaries', 'Angular separation, random distance, 170 binaries', 'Angular separation, random distance, 136 binaries', 
          'Angular separation, distance=100pc, 251 binaries', 'Angular separation, distance=100pc, 170 binaries', 'Angular separation, distance=100pc, 136 binaries']
data = np.roll(data, 6, axis=0)
labels = np.roll(labels, 6, axis=0)
plt.boxplot(np.transpose(data), vert=False, whis=(0, 100), labels=labels)
ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(25))
plt.xticks(np.arange(0.0,1100.0,step=100.0))
plt.xlabel(r'Best fit perturber mass, $M_\odot$')
ax.xaxis.grid(True, which='both', linestyle=':')
plt.xlim([0,1000.0])
plt.show()