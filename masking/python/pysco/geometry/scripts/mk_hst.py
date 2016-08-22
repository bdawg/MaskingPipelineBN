#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt

nr     = 50               # rings within the pupil (should be ~> 50)
rmax   = 2.4/2.           # HST iameter: 2.4 m
rmin   = rmax * 0.33      # central obstruction: (33 %, from TinyTim)
thick  = 0.2              # adopted spider thickness (meters)
srad = 0.15                # square segment size (meters)
rad = np.sqrt(3)*srad     # radius of the first hex ring in meters

xs = np.array(())
ys = np.array(())

fig = plt.figure(0, figsize=(6,6))
plt.clf()
ax = plt.subplot(111)
circ1 = plt.Circle((0,0), rmax, facecolor='none', linewidth=1)
circ2 = plt.Circle((0,0), rmin, facecolor='none', linewidth=1)
ax.add_patch(circ1)
ax.add_patch(circ2)
ax.axis([-rmax,rmax, -rmax,rmax], aspect='equal')

hex_flag = False

if hex_flag:
    for i in range(1-nr, nr, 1):
        for j in xrange(1-nr, nr, 1):
            x = srad * (i + 0.5 * j)
            y = j * np.sqrt(3)/2.*srad
            if (abs(i+j) < nr):
                xs = np.append(xs, x)
                ys = np.append(ys, y)

else:
    for i in range(1-nr, nr, 1):
        for j in xrange(1-nr, nr, 1):
            x = srad * i
            y = srad * j
            xs = np.append(xs, x)
            ys = np.append(ys, y)

# modifications to match the actual telescope pupil (1): diameter constraints
# -----------------------------------------------------------------------
xx, yy = xs.copy(), ys.copy()        # temporary copies
xs, ys = np.array(()), np.array(())  # start from scratch again

osize = 0.05 *rmax # oversize is a fixed % of total telescope size

for i in range(xx.size):
    thisrad = np.sqrt(xx[i]**2 + yy[i]**2)
    if ((rmin+osize) < thisrad < (rmax-osize)):# + 0.1*srad)):
        xs = np.append(xs, xx[i])
        ys = np.append(ys, yy[i])


# modifications to match the actual telescope pupil (2): spiders
# -----------------------------------------------------------
rm_spiders = True

if rm_spiders:
    xx, yy = xs.copy(), ys.copy()        # temporary copies
    xs, ys = np.array(()), np.array(())  # start from scratch again

    epsi = (thick + osize) / np.sqrt(2)
    for i in range(xx.size):
        if ((np.abs(yy[i]) < xx[i]  - epsi) or 
            (np.abs(xx[i]) < yy[i]  - epsi) or
            (xx[i] < -np.abs(yy[i]) - epsi) or
            (yy[i] < -np.abs(xx[i]) - epsi)):
            #(np.abs(yy[i]) > xx[i] + thick/np.sqrt(2))):
        #if ((np.abs(np.arctan(yy[i]/xx[i]) - np.pi/4) > 0.05) and 
        #    (np.abs(np.arctan(yy[i]/xx[i]) + np.pi/4) > 0.05)):
            xs = np.append(xs, xx[i])
            ys = np.append(ys, yy[i])        

ns = xs.size # number of samples


# plot the segments
# ------------------
ax.plot(xs, ys, 'r.')

np.savetxt("../hst2.txt", np.transpose((xs, ys)), 
           fmt='%6.3f')

print "--------------------------------------------------"
print "%d pupil sample points were included in the pupil " % xs.size
print "--------------------------------------------------"
