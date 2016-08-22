''' --------------------------------------------------------------------
  interfere.py - toy routine to read in oifits format data and put this 
  through the kernel phase pipeline

  The oifits.py functions create classes with hidden variables. These 
  are hard to manipulate in pysco and require unwrapping into more 
  manageable arrays.

  The format of giving two baseline coordinates means you need to
  create an inventory of baselines, identify which baselines are 
  referred to by the oifits table, and then build up a transfer matrix
  to map Fourier phases to closure phases.


-------------------------------------------------------------------- '''

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf
import copy
import pickle
import os
import sys
import pdb
import glob
import gzip

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

from scipy.interpolate import griddata

import oifits
import oitools
import emcee
import time
import pysco
import pymultinest

def cp_model(params,u,v,wavel):
    '''Function to model closure phases. Takes a parameter list, u,v triangles and a single wavelength.'''
    ndata = u.shape[0]
    phases = pysco.phase_binary(u.ravel(),v.ravel(),wavel,params)
    phases = np.reshape(phases,(ndata,3))
    cps = np.sum(phases,axis=1)
    return cps

def loglikelihood(params,u,v,wavel,t3data,t3err):
    cps = cp_model(params,u,v,wavel)
    chi2 = np.sum(((t3data-cps)/t3err)**2)
    loglike = -chi2/2
    return loglike


ddir = '/home/ben/Research/Desert/'

oifitsfile = 'Binary1.oifits' #'ScoCen2_Targ1.oifits'

data = oifits.open(ddir+oifitsfile)

ndata = len(data.t3)
for j in data.wavelength:
    wavel = data.wavelength[j].eff_wave
    wavel = wavel[0]
    break

target = data.target[0].target

t3data = []
t3err = []
u = np.zeros((ndata,3))
v = np.zeros((ndata,3))

for j, t3 in enumerate(data.t3):
    t3data.append(t3.t3phi[0])
    t3err.append(t3.t3phierr[0])
    u[j,:] = [t3.u1coord,t3.u2coord,-(t3.u1coord+t3.u2coord)]
    v[j,:] = [t3.v1coord,t3.v2coord,-(t3.v1coord+t3.v2coord)]

t3data = np.array(t3data)
t3err = np.array(t3err)

#----------------------------
# create fake data to test
#----------------------------

mockparams = [90.,220.,8.]
mockerr = t3err*2.
mockdata = cp_model(mockparams,u,v,wavel) + mockerr*np.random.rand(ndata)

'''Now try using the MCMC Hammer for fitting!'''

ivar = np.array([52., 192., 1.53])    # initial parameters for model-fit

ndim, nwalkers = 3, 50

p0 = [ivar + 0.1*ivar*np.random.rand(ndim) for i in range(nwalkers)]

print 'Running emcee now!'

t0 = time.time()

sampler = emcee.EnsembleSampler(nwalkers, ndim, loglikelihood, args=[u,v,wavel,t3data,t3err])
sampler.run_mcmc(p0, 1000)

tf = time.time()

print 'Time elapsed =', tf-t0,'s'

seps = sampler.flatchain[:,0]
ths = sampler.flatchain[:,1]
cs = sampler.flatchain[:,2]

meansep = np.mean(seps)
dsep = np.std(seps)

meanth = np.mean(ths)
dth = np.std(ths)

meanc = np.mean(cs)
dc = np.std(cs)

print 'Separation',meansep,'pm',dsep,'mas'
print 'Position angle',meanth,'pm',dth,'deg'
print 'Contrast',meanc,'pm',dc

paramnames = ['Separation','Position Angle','Contrast']
paramdims = ['(mas)', '(deg)','Ratio']

for i in range(ndim):
    plt.figure(i)
    plt.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
    plt.title(paramnames[i])
    plt.ylabel('Counts')
    plt.xlabel(paramnames[i]+paramdims[i])

plt.show()