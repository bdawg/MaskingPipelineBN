#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt
import pylab
import pyker
from pyker_tools import *

#ddir = './data/GO10879/' # data directory
ddir = './data/GO10143/'

''' ------------------------------------------------------------------

2-frame analysis script.

We want to include both bands of HST data for both the chi-squared and
correlation plots. This script uses pyker instances to manipulate both
frames in parallel.

------------------------------------------------------------------ '''

# --------------------
# 1. load the datasets
# --------------------

num = 25 # index of observation

if ddir == './data/GO10143/':
    Hname = 'n8yj'+str(num)+'010_mos.fits' # data at 1.71 microns
    Jname = 'n8yj'+str(num)+'020_mos.fits' # data at 1.21 microns

elif ddir == './data/GO10879/':
    Hname = 'n9nk'+str(num)+'010_mos.fits'
    Jname = 'n9nk'+str(num)+'020_mos.fits'
else: print 'Data directory failure'

a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "H Band Data"
a.extract_from_fits_frame(ddir + Hname)

b = pyker.KernelPhaseData("./kerphi/hst.pick")
b.name = "J Band Data"
b.extract_from_fits_frame(ddir + Jname)

# --------------------
# 2. load calibrations
# --------------------

a.load_cals('calibs.pick')
b.load_cals('calibs.pick')

if ddir == './data/GO10143/':
    a.calibrate(a.kmean10143H,a.krms10143H)
    b.calibrate(b.kmean10143J,b.krms10143J)
elif ddir == './data/GO10879/':
    a.calibrate(a.kmean10879H,a.krms10879H)
    b.calibrate(b.kmean10879J,b.krms10879J)
else: print 'Data directory failure'

'''--------------------------------------------------------------
Ok, here's the tricky bit. First we guess a pair of contrasts
c1 and c2 and fit over separation to find the position angle
and separation. Then we fix angle and separation and find
c2 and c3 - and fix their ratio. Then we fix angle as well,
and try and tease out the separation-contrast relation.

You may want to pre-optimise the datasets separately! Otherwise
you'll take forever!
--------------------------------------------------------------'''

c1, c2 = 1.49, 2.59

dof = float(a.nkphi + b.nkphi - 3)

# first fix c1, c2 and do sep-angle map

chi2a, xxa, yya = a.fit_map(contrast=c1,sepmin=30,sepmax=100,returnmap=True,
                            display = False)
chi2b, xxb, yyb = b.fit_map(contrast=c2,sepmin=30,sepmax=100,returnmap=True,
                            display = False)
del(xxb)
del(yyb)
chi2 = np.add(chi2a,chi2b)/dof

# find best fits

bestindices = np.unravel_index(chi2.argmin(),chi2.shape)

bestx, besty = xxa[bestindices], yya[bestindices]

sep, angle = besty, bestx

bestfit = chi2.min()

print 'Best chi-squared =', bestfit
print 'Best fit coordinates =', bestx, besty

# --------------------
# optional - plot this
# --------------------

plot = False

if plot == True:
    try:
        clevs = chi2.min()*np.linspace(1,12,20)

        plt.contourf(xxa,yya,chi2,levels = clevs,cmap=plt.cm.bone,
                 antialiased = True)
        plt.axis('tight')
        plt.xlabel('Position Angle (degrees)')
        plt.ylabel('Separation (mas)')
        plt.title('2-Colour Fit Map with c1 ='+str(c1)+', c2 ='+str(c2))
        plt.colorbar()
        plt.draw()
    except: print "/n"

# ---------------------
# now do contrast ratio
# ---------------------

proceed = True

if proceed == True:
    print "Proceeding to contrast map"
    c1min, c1max, c2min, c2max = 1.5,7.0,1.1,10.0
    nx,ny = 120,120

    xaxis = np.linspace(c1min,c1max,num=nx)
    yaxis = np.linspace(c2min,c2max,num=ny)

    xx, yy = np.meshgrid(xaxis,yaxis)

    modelsa = [[phase_binary(a.uv[:,0],a.uv[:,1],a.info['filter'],[sep,angle+90.0 -a.info['orient'],xx[j][i]])
               for i in range(0,nx)] for j in range(0,ny)]

    kphasemodelsa = [[np.dot(a.KerPhi,a.RED*modelsa[j][i]) for i in range (0,nx)]
                    for j in range(0,ny)]

    diffa = [[np.divide(np.subtract(a.kp_signal,kphasemodelsa[j][i])**2,
                       a.kp_error**2) for i in range(0,nx)] for j in range(0,ny)]


    chi2a = [[np.sum(diffa[j][i]) for i in range(0,nx)] for j in range(0,ny)]

    modelsb = [[phase_binary(b.uv[:,0],b.uv[:,1],b.info['filter'],[sep,angle+90.0 -b.info['orient'],yy[j][i]])
               for i in range(0,nx)] for j in range(0,ny)]

    kphasemodelsb = [[np.dot(b.KerPhi,b.RED*modelsb[j][i]) for i in range (0,nx)]
                    for j in range(0,ny)]

    diffb = [[np.divide(np.subtract(b.kp_signal,kphasemodelsb[j][i])**2,
                       b.kp_error**2) for i in range(0,nx)] for j in range(0,ny)]


    chi2b = [[np.sum(diffb[j][i]) for i in range(0,nx)] for j in range(0,ny)]

    chi2 = np.add(chi2a,chi2b)

    chi2 = np.array(chi2)/dof

    bestindices = np.unravel_index(chi2.argmin(),chi2.shape)

    bestx, besty = xx[bestindices], yy[bestindices]

    c1, c2 = bestx, besty

    bestfit = chi2.min()

    print 'Best chi-squared =', bestfit
    print 'Best fit coordinates =', bestx, besty

plot = True

if plot == True:
    try:
        clevs = chi2.min()*np.linspace(1,12,20)

        plt.contourf(xx,yy,chi2,levels = clevs,cmap=plt.cm.bone,
                 antialiased = True)
        plt.axis('tight')
        plt.xlabel('Contrast at '+str(round(a.info['filter']*1e6,2))+' microns')
        plt.ylabel('Contrast at '+str(round(b.info['filter']*1e6,2))+' microns')
        plt.title('2-Colour Fit Map with separation = '+str(round(sep,2))+', position angle = '+str(round(angle,2)))
        plt.colorbar()
        plt.draw()
        plt.show()
    except: print "Failed to plot"
