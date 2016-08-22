#!/usr/bin/env python

''' -------------------------------------------------------
    This procedure generates a coordinates file for a hex
    pupil made of an arbitrary number of rings.
    Additional constraints on the location of spiders make
    it look like the Subaru Telescope primary mirror
    ------------------------------------------------------- '''

import numpy as np, matplotlib.pyplot as plt, pdb
import time

rmax   = 7.92/2.          # outer diameter:      7.92 m
rmin   = 2.30/2.          # central obstruction: 2.30 m
thick  = 0.25             # adopted spider thickness (meters)
offset = 1.278            # spider intersection offset (meters)
beta   = 51.75*np.pi/180  # spider angle beta
epsi   = thick/(2*np.sin(beta))

nr     = 50               # number of rings within the pupil

srad = 0.39            # segment "radius" (diagonal = 1.8 m)
rad = np.sqrt(3)*srad # radius of the first hex ring in meters

no_spider = False#True#


xs = np.array(())
ys = np.array(())

fig = plt.figure(figsize=(6,6))
ax = plt.subplot(111)
circ1 = plt.Circle((0,0), rmax, facecolor='none', linewidth=1)
circ2 = plt.Circle((0,0), rmin, facecolor='none', linewidth=1)
ax.add_patch(circ1)
ax.add_patch(circ2)
#plt.clf()
ax.axis([-rmax,rmax, -rmax,rmax], aspect='equal')

for i in range(1-nr, nr, 1):
    for j in xrange(1-nr, nr, 1):
        x = srad * (i + 0.5 * j)
        y = j * np.sqrt(3)/2.*srad
        if (abs(i+j) < nr):
            xs = np.append(xs, x)
            ys = np.append(ys, y)
            #plt.text(x, y, '%d,%d' % (i,j))    

# modifications to match the actual Subaru pupil (1): diameter constraints
# -----------------------------------------------------------------------
xx, yy = xs.copy(), ys.copy()        # temporary copies
xs, ys = np.array(()), np.array(())  # start from scratch again

for i in range(xx.size):
    thisrad = np.sqrt(xx[i]**2 + yy[i]**2)
    if (1.05 * rmin < thisrad < (0.95 * rmax)):# + 0.1*srad)):
        xs = np.append(xs, xx[i])
        ys = np.append(ys, yy[i])

# plot the freakin spiders
# ------------------------
#plt.plot([offset/2.,rmax], [0.0, rmax+offset/2*np.cos(beta)], 'b', lw=5)
#plt.plot([offset/2.,rmax], [0.0, -rmax-offset/2*np.cos(beta)], 'b', lw=5)
#plt.plot([-offset/2.,-rmax], [0.0, rmax+offset/2*np.cos(beta)], 'b', lw=5)
#plt.plot([-offset/2.,-rmax], [0.0, -rmax-offset/2*np.cos(beta)], 'b', lw=5)

# modifications to match the actual Subaru pupil (2): spiders
# -----------------------------------------------------------
rm_spiders = False

if rm_spiders:
    xx, yy = xs.copy(), ys.copy()        # temporary copies
    xs, ys = np.array(()), np.array(())  # start from scratch again

    for i in range(xx.size):
        xt1, xt2 =  xx[i] - offset/2., xx[i] + offset/2.
        th1 = np.arctan(yy[i]/(xt1 - epsi)) 
        th2 = np.arctan(yy[i]/(xt2 + epsi))
        th3 = np.arctan(yy[i]/(xt1 + epsi)) 
        th4 = np.arctan(yy[i]/(xt2 - epsi)) 

        if ((xt1 > epsi) and (np.abs(th1) < beta)):
            xs = np.append(xs, xx[i])
            ys = np.append(ys, yy[i])        

        if ((xt2 < -epsi) and (np.abs(th2) < beta)):
            xs = np.append(xs, xx[i])
            ys = np.append(ys, yy[i])        

        if ((yy[i] > 0) and ((np.abs(th3) > beta) or np.abs(th4) > beta)):
            xs = np.append(xs, xx[i])
            ys = np.append(ys, yy[i])        
        
        if ((yy[i] < 0) and ((np.abs(th3) > beta) or np.abs(th4) > beta)):
            xs = np.append(xs, xx[i])
            ys = np.append(ys, yy[i])        

# plot segments
# -------------
r0 = srad/np.sqrt(3)
th = 2*np.pi*np.arange(6)/6. + np.pi/6.

for i in range(xs.size):
    hx = xs[i] + r0 * np.cos(th)
    hy = ys[i] + r0 * np.sin(th)
    ax.fill(hx, hy, fc='none', linewidth=1)

ax.plot(xs, ys, 'r.')
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

#np.savetxt("../subaru_hex_342.txt", np.transpose((xs,ys)), 
#           fmt='%12.9f')


np.savetxt("../test.txt", 
           np.transpose((np.delete(xs,np.where((ys == 0) * (xs <0))),
                         np.delete(ys,np.where((ys == 0) * (xs <0))))),
           fmt='%12.9f')

print "--------------------------------------------------"
print "%d pupil sample points were included in the pupil " % xs.size
print "--------------------------------------------------"



plt.show()

