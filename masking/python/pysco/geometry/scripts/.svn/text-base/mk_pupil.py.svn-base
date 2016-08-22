#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt

nb_sub = 5  # number of sub-apertures accross one pupil-diameter
obstr  = 0.0 # size of central obstruction

x = np.linspace(0, 2, nb_sub) - 1.0
y = np.linspace(0, 2, nb_sub) - 1.0

coords = np.zeros((0,2))

for i in range(nb_sub):
    for j in range(nb_sub):
        #if obstr <= np.hypot(x[i],y[j]) <= 1.0:
        coords=np.append(coords, [[x[i], y[j]]], axis=0)

plt.clf()
plt.plot(coords[:,0], coords[:,1], 'ro')
plt.axis([-2,2,-2,2])
np.savetxt("coords.txt", np.transpose((coords[:,0], coords[:,1])), 
           fmt='%6.3f')
plt.show()
