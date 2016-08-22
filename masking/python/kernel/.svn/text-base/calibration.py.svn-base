#!/usr/bin/env python

import pyker
import numpy as np
import matplotlib.pyplot as plt

#ddir = './data/GO10879/' # data directory
ddir = './data/GO10143/'

''' ------------------------------------------------------------------
This module is meant to calibrate kernel phases in HST data.

We select those sources we believe to be pointlike, and then
calculate kernel phases in each case. We then find the mean and
uncertainty on each kernel phase and use this as our calibration.
----------------------------------------------------------------- '''

#create lists of filenames

GO10143_1 = ['n8yj'+str(j)+'010_mos.fits' for j in ['0'+str(k) for k in range(1,10)] + range(10,69)] #GO10143 obs at 1.71 um
GO10143_2 = ['n8yj'+str(j)+'020_mos.fits' for j in ['0'+str(k) for k in range(1,10)] + range(10,69)] #GO10143 at 1.12 um

GO10879_1 = ['n9nk'+str(j)+'010_mos.fits' for j in ['0'+str(k) for k in range(1,10)] + range(10,32)] #GO10879 at 1.71 um
GO10879_2 = ['n9nk'+str(j)+'020_mos.fits' for j in ['0'+str(k) for k in range(1,10)] + range(10,32)] #GO10879 at 1.12 um

# -------------------
# 1. load the dataset
# -------------------

data = GO10143_1
#data = GO10143_2
#data = GO10879_1
#data = GO10879_2

#data = [GO10143_1,GO10143_2,GO10879_1,GO10879_2] # for looping over datasets

''' ------------------------------------------------------------------
Calibration sources are selected from my table. These are lists of sources
for which no correlation or visual abnormality was noted.

CAUTION: Absence of detection is != detection of absence!
----------------------------------------------------------------- '''

caldata = [27,31,32,33,35,36,42,45,49,51,52,53,54,56,57,60,63,65,66,67] #GO10143
#caldata = [1,6,8,11,13,15,17,18,19,20,22,23,25,26,29,31] # GO10879

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst20.pick")
a.name = "HST - 10879 trial" # add a label to the template
num = a.nkphi

# load the FITS frames, and extract the Kernel-phases using the already
# loaded HST template.

kerphases = np.zeros((len(caldata),a.nkphi))

for n,j in enumerate(caldata):
    print ddir + data[j-1]
    a.extract_from_fits_frame(ddir + data[j-1])
    kerphases[n,:] = a.kp_signal

#calculate mean, mean square difference from mean, rms

kmean = np.array([np.mean([(kerphases[j-1])[k] for j in range(0,len(caldata))]) for k in range(0,num)])

kmsd= np.array([np.mean((([(kerphases[j-1])[k] for j in range(0,len(caldata))])-kmean[k])**2) for k in range(0,num)])

krms = np.sqrt(kmsd)


#now plot histograms
##
##plt.clf()
##
##f0 = plt.subplot(121)
##f0.hist(kmean,30)
##plt.xlabel('Kernel Phase Ensemble Mean (deg)')
##plt.ylabel('Occurrences')
##plt.title('Calibration GO10143 2')
##f1 = plt.subplot(122)
##f1.hist(krms,30)
##plt.xlabel('Kernel Phase Ensemble RMS error (deg)')
##plt.ylabel('Occurrences')
##plt.title('RMS Error')
##plt.draw()
##plt.savefig('hist143_2')
##plt.show()
