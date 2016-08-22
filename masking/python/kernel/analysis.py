#!/usr/bin/env python

import pyker
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyker_tools import *

#ddir = './data/GO10879/' # data directory
ddir = "./data/GO10143/"

pickdir = './chi2/'

num = 11 # index of file in database
band = 1 # H = 1, J = 2



''' ------------------------------------------------------------------
Basic 1-frame analysis script.

After the "import pyker" command, the documentation for the code is
available typing "help(pyker)".

Note that in the calibrations we refer to dataset 2 as J, as it consists
of J band observations, and to dataset 1 as H likewise.

There are indexed lists of filenames below in case you want to perform
analysis on a batch of files using a loop; otherwise it's easy to just
change the name of the file. Careful to do the directory.

----------------------------------------------------------------- '''
#create lists of filenames
#GO10143_1 = ['n8yj'+str(j)+'010_mos.fits' for j in ['0'+str(k) for k in range(1,9)] + range(0,68)] #GO10143 obs at 1.71 um
#GO10143_2 = ['n8yj'+str(j)+'020_mos.fits' for j in ['0'+str(k) for k in range(1,9)] + range(0,68)] #GO10143 at 1.12 um

#GO10879_1 = ['n9nk'+str(j)+'010_mos.fits' for j in ['0'+str(k) for k in range(1,9)] + range(10,68)] #GO10873 at 1.71 um
#GO10879_2 = ['n9nk'+str(j)+'020_mos.fits' for j in ['0'+str(k) for k in range(1,9)] + range(10,68)] #GO10873 at 1.12 um

# -------------------
# 1. load the dataset
# -------------------

if band==1:
    bandn = 'H'
elif band==2:
    bandn = 'J'
else: print 'Invalid band'

if num >= 10:
    num = str(num)
else:
    num = '0'+str(num)
    
if ddir == './data/GO10143/':
    fname = 'n8yj'+num+'0'+str(band)+'0_mos.fits'
elif ddir == './data/GO10879/':
    fname = 'n9nk'+num+'0'+str(band)+'0_mos.fits'
else: print 'Data directory failure'

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "HST - " + str(num)+bandn #  # add a label to the template

# load the FITS frame, and extract the Kernel-phases using the already
# loaded HST template.

#a.extract_from_fits_frame(ddir + 'n9nk28010_mos.fits') # for selecting filenames manually
#a.extract_from_fits_frame(ddir + 'n8yj25010_mos.fits')

a.extract_from_fits_frame(ddir + fname)

#a.extract_from_fits_frame(ddir + GO10143_1[5]) # for convenience
#a.extract_from_fits_frame(ddir + GO10879_1[5])

# --------------------------------------------------------
# a couple of optional plots which you may want to look at
# --------------------------------------------------------
#a.plot_pupil_and_uv(1.5) # optional plot of Kernel-phase model
#a.plot_im_and_ac()       # optional plot of image and Fourier counterpart

# --------------------
# 2. load calibrations
# --------------------

#a.calibrate(0,np.std(a.kp_signal)) 

a.load_cals('calibs.pick')
if ddir == './data/GO10143/':
    if band == 1:
        a.calibrate(a.kmean10143H,a.krms10143H)
    elif band == 2:
        a.calibrate(a.kmean10143J,a.krms10143J)
    else: print 'Band failure'
elif ddir == './data/GO10879/':
    if band == 1:
        a.calibrate(0,a.krms10879H)
    elif band == 2:
        a.calibrate(0,a.krms10879J)
    else: print 'Band failure'
else: print 'Data directory failure'

# ------------------------
# 2. model fit of the data
# ------------------------

# binary_fit uses Levenberg-Marquardt Algorithm (LMA) to minimize the
# variance between the data and a binary-star model. Just like every
# implementation of LMA, you need to provide a initial set of
# parameters.
# it returns the best fit parameters, and a covariance matrix that can
# be used to determine uncertainties on the fit.

##print "up to model fitting"
##init_params = [70, 128,0.88]    # initial parameters for model-fit
##optim = a.binary_fit(init_params)  # fit data with binary model
##params = optim[0]                  # best fit parameters
##
##b = pyker.KernelPhaseData("./kerphi/hst.pick")
##b.binary_model([48,128,0.88], a.info)
##
# -------------------
# 3. correlation plot
# -------------------

#a.correlation_plot(b)
#a.residual_plot(b)

##if optim[1] != None:
##    print "\nParameter estimate covariance matrix:"
##    print "-----------------------------------\n"
##    print np.round((a.rms)**2 * optim[1], 4)
##
##if optim[0] != None:
##    print "\nParameter estimate:"
##    print "-----------------------------------\n"
##    print optim[0]
##plt.show()
#plt.savefig('cplt25')

#a.fit_map(angle=120,sepmin=30,sepmax=100,cmin=3.0,cmax=70,nclevs=40)

##sepmin,sepmax,anglemin,anglemax,cmin,cmax = 40.,90.,115.,150.,1.1,2.0
##nsep,nangle,nc = 200,200,300
##
##save_chi2(a,pickdir+str(num)+bandn+'.pick',sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc)
##
##gridH = fully_marginalise(*load_chi2('21H.pick'))
##gridJ = fully_marginalise(*load_chi2('21J.pick'))
##
##plotnum = 0
##
##plt.clf()
##plt.plot(gridH[plotnum],gridH[plotnum+3]/max(gridH[plotnum+3]),
##         gridJ[plotnum],gridJ[plotnum+3]/max(gridJ[plotnum+3]))
##plt.ylabel('Likelihood')
##plt.title('Marginalised Likelihood Function')
##plt.legend(('H Band','J Band'))
##
##if plotnum == 0: plt.xlabel('Separation (mas)')
##elif plotnum == 1: plt.xlabel('Position Angle (degrees')
##elif plotnum == 2: plt.xlabel('Contrast Ratio')
##else: print 'Incorrect plot number'
##
##plt.show()
##
params = np.array([70.,50.,10.0])

ranges = np.array([0.1,1,0.07])*1.5

distribution = a.mh(params,ranges,1e6,cmax=80)
