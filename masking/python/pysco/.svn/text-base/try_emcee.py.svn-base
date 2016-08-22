#!/usr/bin/env python

import pysco
import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
import pysco.fitting as fit
import emcee
import time

# provide an absolute path to the data directory
# ----------------------------------------------
home = os.environ['HOME']
ddir = './data/' # data directory

''' ------------------------------------------------------------------
A hack of Frantz' try_nicmos.py to test out the emcee 'MCMC Hammer' 
algorithm. 
------------------------------------------------------------------ '''

def lnprob(x, kpo1,kpo2):
    x1 = [x[0],x[1],x[2]]
    x2 = [x[0],x[1],x[3]]

    modl_ker1 = fit.binary_model(x1,kpo1.kpi,kpo1.hdr)
    modl_ker2 = fit.binary_model(x2,kpo2.kpi,kpo2.hdr)

    chisquared = np.sum(((modl_ker1 - kpo1.kpd)/(kpo1.kpe))**2)+np.sum(((modl_ker2 - kpo2.kpd)/(kpo2.kpe))**2)
    
    return -chisquared/2.


# -------------------------------
# 1. create the KP info structure
# -------------------------------

# once saved, the kpi.gz structure can be directly reloaded when 
# creating a KPO instance, such as done in step #2.

a = pysco.KPI("./geometry/nicmos.txt")
a.name = "HST - NIC1" #  # add a label to the template
a.save_to_file('./hst.kpi.gz')

# -------------------
# 2. load the dataset
# -------------------

# load the FITS frame, and extract the Kernel-phases using the
# HST KPI template calculated at the previous step. 
# Two data sets are provided:
# n8yj59010_mos.fits.gz and 'n8yj59020_mos.fits.gz

a = pysco.KPO('./hst.kpi.gz')
b = pysco.KPO('./hst.kpi.gz')
a.extract_KPD(ddir+'n8yj59010_mos.fits', plotim=False)
b.extract_KPD(ddir+'n8yj59020_mos.fits', plotim=False)

# ---------------------------------------

# ------------------------
# 2. model fit of the data
# ------------------------

ivar = np.array([150.0, 80.0, 5.0,5.0])    # initial parameters for model-fit

ndim, nwalkers = 4, 100

p0 = [ivar + np.random.rand(ndim) for i in range(nwalkers)]

print 'Running emcee now!'

t0 = time.time()

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[a,b])
sampler.run_mcmc(p0, 1000)

tf = time.time()

print 'Time elapsed =', tf-t0,'s'

seps = sampler.flatchain[:,0]
ths = sampler.flatchain[:,1]
c1s = sampler.flatchain[:,2]
c2s = sampler.flatchain[:,3]

meansep = np.mean(seps)
dsep = np.std(seps)

meanth = np.mean(ths)
dth = np.std(ths)

meanc1 = np.mean(c1s)
dc1 = np.std(c1s)

meanc2 = np.mean(c2s)
dc2 = np.std(c2s)

print 'Separation',meansep,'pm',dsep,'mas'
print 'Position angle',dth,'pm',dth,'deg'
print 'H Contrast',meanc1,'pm',dc1
print 'J contrast',meanc2,'pm',dc2

# for i in range(ndim):
#     plt.figure(i)
#     plt.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
#     plt.title("Dimension {0:d}".format(i))

# plt.show()