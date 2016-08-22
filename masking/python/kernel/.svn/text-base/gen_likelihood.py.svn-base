#!/usr/bin/env python

import pyker
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyker_tools import *
plt.ioff()
ddir43 = './data/GO10143/'
ddir79 = './data/GO10879/'

targets = ['25H','25J','21H','21J','47H']

fnames = [ddir43+'n8yj25010_mos.fits',ddir43+'n8yj25020_mos.fits',
         ddir79+'n9nk21010_mos.fits',ddir79+'n9nk21020_mos.fits',
          ddir43+'n8yj47010_mos.fits']

chi2s = ['25Hdummy.pick','25Jlower.pick','21H.pick','21J.pick','47H.pick']

for j,name in enumerate(fnames):
    print j,name    
    a = pyker.KernelPhaseData("./kerphi/hst.pick")
    a.name = targets[j]
    a.extract_from_fits_frame(name)
    chi2,sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc = load_chi2(chi2s[j])

    for i in [0,1,2]:
        print i
        plt.clf()
        a.marginalplot(i,sepmin,sepmax,anglemin,anglemax,cmin,cmax,
                       nsep=nsep,nangle=nangle,nc=nc,chi2=chi2,
                       save='./Likelihood/'+targets[j]+str(i)+'.png')

    
