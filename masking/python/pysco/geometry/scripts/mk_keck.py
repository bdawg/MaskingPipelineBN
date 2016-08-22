#!/usr/bin/env python

''' -------------------------------------------------------
    This procedure generates a coordinates file for Ker-phi
    analysis of data acquired at Keck. Four steps:
    - coordinates of the center of each segment
    - generate a hex-pattern per segment
    - eliminate all doublons
    - remove coordinates along spider arms
    ------------------------------------------------------- '''

import numpy as np, matplotlib.pyplot as plt, pdb

srad = 0.9 # segment "radius" (diagonal = 1.8 m)
rad = np.sqrt(3)*srad # radius of the first hex ring in meters

no_spider = False#True#

nr = 6

th = 2*np.pi*np.arange(nr)/nr

# -----------------------------------------------------------------------
# here, we generate a series of coordinates that correspond to the center
# of each segment.
# -----------------------------------------------------------------------
# inner ring
x = rad*np.cos(th + np.pi/6.)
y = rad*np.sin(th + np.pi/6.)

# second ring
x = np.append(x, 2*rad*np.cos(th+np.pi/6.))
y = np.append(y, 2*rad*np.sin(th+np.pi/6.))
x = np.append(x, np.sqrt(3)*rad*np.cos(th))
y = np.append(y, np.sqrt(3)*rad*np.sin(th))

#outer ring
x = np.append(x, 3*rad*np.cos(th+np.pi/6.))
y = np.append(y, 3*rad*np.sin(th+np.pi/6.))

r0 = np.sqrt(21)*srad
x = np.append(x, r0*np.cos(th+np.arctan(2./np.sqrt(3))))
y = np.append(y, r0*np.sin(th+np.arctan(2./np.sqrt(3))))

x = np.append(x, r0*np.cos(th+np.arctan(2./np.sqrt(3))))
y = np.append(y, -r0*np.sin(th+np.arctan(2./np.sqrt(3))))

# -----------------------------------------------------------------------
# now these coordinates can be used to generate the preliminary sampling
# map: coordinates + off-center (hex pattern also) to locate samples
# -----------------------------------------------------------------------

srad = 0.45 # sample point separation along a "diagonal"
rad = np.sqrt(3)*srad # radius of the first hex ring in meters

xx, yy = 0, 0

for i in range(x.size):
    # "inner ring"
    addx = np.append(x[i] + srad*np.cos(th), x[i])
    addy = np.append(y[i] + srad*np.sin(th), y[i])
    # "second ring"
    addx = np.append(addx, x[i] + 2*srad*np.cos(th))
    addy = np.append(addy, y[i] + 2*srad*np.sin(th))

    addx = np.append(addx, x[i] + rad*np.cos(th+np.pi/6.))
    addy = np.append(addy, y[i] + rad*np.sin(th+np.pi/6.))

    xx = np.append(xx, addx)
    yy = np.append(yy, addy)

xx = np.delete(xx, 0)
yy = np.delete(yy, 0)

xx = x
yy = y

# -----------------------------------------------------------------------
# need to elimitate all double entries in the coordinates introduced by
# the previous procedures
# -----------------------------------------------------------------------
ndgt = 2 # number of digits of precision for rounding
prec = 10**(-ndgt)

a = np.unique(np.round(xx, ndgt)) # distinct u-component of baselines
nbx = a.shape[0]                   # number of distinct u-components

xx_sel, yy_sel = 0,0 # empty array to store selected baselines

for i in range(nbx):     # figure out distinct v-components and fill uv_sel
    b = np.where(np.abs(xx - a[i]) <= prec)
    c = np.unique(np.round(yy[b], ndgt))
    nby = np.shape(c)[0] # number of distinct v-compoments
    for j in range(nby):
        xx_sel = np.append(xx_sel, a[i])
        yy_sel = np.append(yy_sel, c[j])

xx = xx_sel
yy = yy_sel

xx = np.delete(xx, 0)
yy = np.delete(yy, 0)


# -----------------------------------------------------------------------
# eliminate sample points that lie on the spider arms, conveniently
# located at azimuths multiples of 60 degres!
# -----------------------------------------------------------------------

if no_spider:
    a = (abs(np.arctan2(yy,xx)) <= 0.01)
    diag = a
    a = (abs(np.arctan2(yy,xx) - np.pi/3.) <= 0.01)
    diag += a
    a = (abs(np.arctan2(yy,xx) - 2*np.pi/3.) <= 0.01)
    diag += a
    a = (abs(np.arctan2(yy,xx) - np.pi) <= 0.01)
    diag += a
    a = (abs(np.arctan2(yy,xx) + 2*np.pi/3.) <= 0.01)
    diag += a
    a = (abs(np.arctan2(yy,xx) + np.pi/3.) <= 0.01)
    diag += a

    keep = 1 - diag

    xx = xx[np.where(keep)]
    yy = yy[np.where(keep)]

# -----------------------------------------------------------------------
# finally, plot the pupil sampling map and save the coordinates in a
# text file for the Ker-phase generation program.
# -----------------------------------------------------------------------

plt.clf()


# plot segments
# -------------
r0 = 0.9
for i in range(x.size):
    hx = x[i] + r0 * np.cos(th)
    hy = y[i] + r0 * np.sin(th)
    plt.fill(hx, hy, fc='none', linewidth=1)

hx = r0 * np.cos(th)
hy = r0 * np.sin(th)

plt.fill(hx, hy, color='gray')

# plot spider arms
# ----------------
r0 = 5 * 0.9
x0, x1, y0, y1 = -r0, r0, 0, 0
plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)
x0, y0 = r0 * np.cos(np.pi/3), r0 * np.sin(np.pi/3)
x1, y1 = -x0, -y0
plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)
x0, y0 = r0 * np.cos(2*np.pi/3), r0 * np.sin(2*np.pi/3)
x1, y1 = -x0, -y0
plt.plot([x0, x1], [y0, y1], color='gray', linewidth=5)


# plot sample points
# ------------------
plt.plot(xx, yy, 'ro')
plt.axis([-6,6, -6,6], aspect='equal')

np.savetxt("geometry/keck1.txt", np.transpose((xx,yy)), 
           fmt='%6.3f')

print "--------------------------------------------------"
print "%d pupil sample points were included in the pupil " % xx.size
print "--------------------------------------------------"

plt.show()
