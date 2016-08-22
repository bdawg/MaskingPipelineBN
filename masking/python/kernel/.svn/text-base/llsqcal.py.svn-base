#!/usr/bin/env python

import pyker
import numpy as np
import pylab as plt
import pickle
from pyker_tools import *

#ddir = './data/GO10879/' # data directory
ddir = './data/GO10143/'

num = 25 # index of file in database
band = 1 # H = 1, J = 2

if band == 1:
    bandname = 'H'
elif band == 2:
    bandname = 'J'
else: print 'Invalid band!'

''' ------------------------------------------------------------------

Toy for manipulating linear least squares calibration.

----------------------------------------------------------------- '''

# unpickle calibration operators

myf = open('calops.pick','r')
calops = pickle.load(myf)
myf.close()

if num >= 10: # ensure the number is padded with zeros in the fname
    num = str(num)
else:
    num = '0'+str(num)
    
if ddir == './data/GO10143/':
    fname = 'n8yj'+num+'0'+str(band)+'0_mos.fits'
    A = calops['op143'+bandname]
    
elif ddir == './data/GO10879/':
    fname = 'n9nk'+num+'0'+str(band)+'0_mos.fits'
    A = calops['op879'+bandname]
else: print 'Data directory failure'

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "HST - 10143 trial" #  # add a label to the template

a.extract_from_fits_frame(ddir + fname)

# calibrate a with linear least squares

systematic, residues, rank, sv = np.linalg.lstsq(A,a.kp_signal)

correction = np.dot(A,systematic)

#print correction
#print correction.shape

# now apply the systematic calibration

a.corr_signal = a.kp_signal - np.dot(A,systematic)

b = pyker.KernelPhaseData("./kerphi/hst.pick")

params = [70.,349.,1.49]

b.binary_model(params,a.info)

sysb,resb,rankb,svb = np.linalg.lstsq(A,b.kp_signal)

correctb = np.dot(A,sysb)
