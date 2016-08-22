#!/usr/bin/env python

import pyker
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyker_tools import *

ddir43 = './data/GO10143/'
ddir79 = './data/GO10879/'

a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "Candidate 21 H band"
a.extract_from_fits_frame(ddir79 + 'n9nk21020_mos.fits')
a.load_cals('calibs.pick')
a.calibrate(a.kmean10879H,a.krms10879H)

sepmin,sepmax,anglemin,anglemax,cmin,cmax = 30.,80.,95.,120.,1.05,2.5
nsep,nangle,nc = 200,200,300

save_chi2(a,'21H.pick',sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc)

##del(a)
##
##a = pyker.KernelPhaseData("./kerphi/hst.pick")
##a.name = "Candidate 47 J band"
##a.extract_from_fits_frame(ddir43 + 'n8yj47020_mos.fits')
##a.load_cals('calibs.pick')
##a.calibrate(a.kmean10143J,a.krms10143J)
##
##sepmin,sepmax,anglemin,anglemax,cmin,cmax = 50.,70.,125.,140.,1.3,3.0
##nsep,nangle,nc = 200,200,200
##
##save_chi2(a,'47J.pick',sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc)
