#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt

nb_sub = 21  # number of sub-apertures accross one pupil-diameter

prad = 1.0 # pupil radius in meters
prad = 0.8 # scaled down 10 %

x = np.linspace(0, 2*prad, nb_sub) - prad
y = np.linspace(0, 2*prad, nb_sub) - prad

xi,yi = np.meshgrid(x, y)
dist  = np.hypot(yi,xi)

hmax = 2.0*prad/3.0*np.exp(-np.abs(2*x/prad)**2.0)

coords = np.zeros((0,2))

for i in range(nb_sub):
    for j in range(nb_sub):
        if np.abs(yi[i,j]) <= hmax[j]:
            coords=np.append(coords, [[x[i], y[j]]], axis=0)
        #if dist[i,j] <= prad:
        #    coords=np.append(coords, [[x[i], y[j]]], axis=0)
plt.clf()
plt.plot(coords[:,0], coords[:,1], 'ro')
plt.axis([-2,2,-2,2])
np.savetxt("geometry/shaped_pup.txt", np.transpose((coords[:,0], coords[:,1])), 
           fmt='%6.3f')
#plt.show()
