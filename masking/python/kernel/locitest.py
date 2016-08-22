#!/usr/bin/env python

import pyker
import numpy as np
import pickle
import matplotlib.pyplot as plt

ddir = './data/GO10879/' # data directory
#ddir = './data/GO10143/'

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

#data = GO10143_1
#data = GO10143_2
data = GO10879_1
#data = GO10879_2

#data = [GO10143_1,GO10143_2,GO10879_1,GO10879_2] # for looping over datasets

''' ------------------------------------------------------------------
Calibration sources are selected from my table. These are lists of sources
for which no correlation or visual abnormality was noted.

CAUTION: Absence of detection is != detection of absence!
----------------------------------------------------------------- '''

#caldata = [26,27,31,32,33,35,36,42,43,45,49,51,52,53,54,56,57,60,63,65,66,67] #GO10143
caldata = [1,6,8,11,13,15,17,18,19,20,22,23,25,26,29,31] # GO10879

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "HST - 10879 trial" # add a label to the template


# load the FITS frames, and extract the Kernel-phases using the already
# loaded HST template.

kerphases879H = np.empty((a.nkphi,len(caldata)))

for j,datanum in enumerate(caldata):
    print ddir + data[datanum-1]
    a.extract_from_fits_frame(ddir + data[datanum-1])
    kerphases879H[:,j] = a.kp_signal

# now do J band

data = GO10879_2

kerphases879J = np.empty((a.nkphi,len(caldata)))

for j,datanum in enumerate(caldata):
    print ddir + data[datanum-1]
    a.extract_from_fits_frame(ddir + data[datanum-1])
    print j
    kerphases879J[:,j] = a.kp_signal

# now do H band for 10143

ddir = './data/GO10143/'
data = GO10143_1

caldata = [26,27,31,32,33,35,36,42,43,45,49,51,52,53,54,56,57,60,63,65,66,67]

kerphases143H = np.empty((a.nkphi,len(caldata)))

for j,datanum in enumerate(caldata):
    print ddir + data[datanum-1]
    a.extract_from_fits_frame(ddir + data[datanum-1])
    kerphases143H[:,j] = a.kp_signal

# now do J band for 10143

data = GO10143_2

kerphases143J = np.empty((a.nkphi,len(caldata)))

for j,datanum in enumerate(caldata):
    print ddir + data[datanum-1]
    a.extract_from_fits_frame(ddir + data[datanum-1])
    kerphases143J[:,j] = a.kp_signal

# now pickle these structures

file = 'calops.pick'

try: cals = {'op143H'   :   kerphases143H,
             'op143J'   :   kerphases143J,
             'op879H'   :   kerphases879H,
             'op879J'   :   kerphases879J}
except:
    print('Calibration failed!')

try: myf = open(file,'w')
except:
    print 'File cannot be created.'

pickle.dump(cals,myf)
myf.close()
