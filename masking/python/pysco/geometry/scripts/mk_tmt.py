#!/usr/bin/env python

''' -------------------------------------------------------
    This procedure generates a coordinates file for a hex
    pupil made of an arbitrary number of rings.
    ------------------------------------------------------- '''

import numpy as np, matplotlib.pyplot as plt, pdb
import time

srad = 1.3 # segment "radius" (diagonal = 1.8 m)
rad = np.sqrt(3)*srad # radius of the first hex ring in meters

no_spider = False#True#

nr = 14
rmax = np.round(nr * srad)

xs = np.array(())
ys = np.array(())

plt.clf()
plt.axis([-rmax,rmax, -rmax,rmax], aspect='equal')

for i in range(1-nr, nr, 1):
    for j in xrange(1-nr, nr, 1):
        x = srad * (i + 0.5 * j)
        y = j * np.sqrt(3)/2.*srad
        if (abs(i+j) < nr):
            xs = np.append(xs, x)
            ys = np.append(ys, y)
            #plt.text(x, y, '%d,%d' % (i,j))    

# modifications to match the actual TMT pupil
xx, yy = xs.copy(), ys.copy()
xs, ys = np.array(()), np.array(())

for i in range(xx.size): 
    if (srad < np.sqrt(xx[i]**2 + yy[i]**2) < rmax*0.97*np.sqrt(3)/2.):
#    if (10.7 * srad < np.sqrt(xx[i]**2 + yy[i]**2) < rmax*0.97*np.sqrt(3)/2.):
        xs = np.append(xs, xx[i])
        ys = np.append(ys, yy[i])

# plot segments
# -------------
r0 = srad/sqrt(3)
th = 2*np.pi*np.arange(6)/6. + np.pi/6.

for i in range(xs.size):
    hx = xs[i] + r0 * np.cos(th)
    hy = ys[i] + r0 * np.sin(th)
    plt.fill(hx, hy, fc='none', linewidth=1)

plt.plot(xs, ys, 'r.')
# plot spider arms
# ----------------
#r0 = rmax
#x0, x1, y0, y1 = -r0, r0, 0, 0
#plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)
#x0, y0 = r0 * np.cos(np.pi/3), r0 * np.sin(np.pi/3)
#x1, y1 = -x0, -y0
#plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)
#x0, y0 = r0 * np.cos(2*np.pi/3), r0 * np.sin(2*np.pi/3)
#x1, y1 = -x0, -y0
#plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)


# plot sample points
# ------------------
#plt.plot(xx, yy, 'ro')
#plt.axis([-6,6, -6,6], aspect='equal')

np.savetxt("../tmt.txt", np.transpose((xs,ys)), 
           fmt='%12.9f')

print "--------------------------------------------------"
print "%d pupil sample points were included in the pupil " % xs.size
print "--------------------------------------------------"

plt.show()
