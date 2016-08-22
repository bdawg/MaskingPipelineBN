#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt

rad = 1.0 # radius of the first hex ring in meters

nr = 6

th = 2*np.pi*np.arange(nr)/nr

# inner ring
x = np.cos(th + np.pi/6.)
y = np.sin(th + np.pi/6.)

# outer ring
x = np.append(x, 2*rad*np.cos(th+np.pi/6.))
y = np.append(y, 2*rad*np.sin(th+np.pi/6.))
x = np.append(x, np.sqrt(3)*rad*np.cos(th))
y = np.append(y, np.sqrt(3)*rad*np.sin(th))


plt.clf()
plt.plot(x, y, 'ro')
plt.axis([-3,3, -3,3])

#np.savetxt("jwst.txt", np.transpose((x,y)), 
#           fmt='%6.3f')
plt.show()
