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

# -------------------------------
# 1. create the KP info structure
# -------------------------------

# once saved, the kpi.gz structure can be directly reloaded when 
# creating a KPO instance, such as done in step #2.

# a = pysco.KPI("./geometry/nicmos.txt")
# a.name = "HST - NIC1" #  # add a label to the template
# a.save_to_file('./hst.kpi.gz')


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
a.kpi.name = "2M XXXX-XX" #  # labels the data
b.kpi.name = '2M XXXX-XX'

# ---------------------------------------

# ------------------------
# 2. model fit of the data
# ------------------------

ivar = np.array([150.0, 80.0, 5.0,5.0])    # initial parameters for model-fit

ndim, nwalkers,ntemps = 4, 100,4

p0 = np.random.uniform(low=-1.0, high=1.0, size=(ntemps, nwalkers, ndim))

for j in range(ntemps):
    for k in range(nwalkers):
        p0[j,k,:] += ivar

def lnprob(x):
    x1 = [x[0],x[1],x[2]]
    x2 = [x[0],x[1],x[3]]

    modl_ker1 = fit.binary_model(x1,a.kpi,a.hdr)
    modl_ker2 = fit.binary_model(x2,b.kpi,b.hdr)

    chisquared = np.sum(((modl_ker1 - a.kpd)/(a.kpe))**2)+np.sum(((modl_ker2 - b.kpd)/(b.kpe))**2)
    
    return -chisquared/2.
    
def lnprior(x):
    return 0.0

print 'Running emcee now!'

t0 = time.time()

sampler = emcee.PTSampler(ntemps,nwalkers, ndim, lnprob,lnprior,threads=4)
for p, lnprob, lnlike in sampler.sample(p0, iterations=1000):
    pass
sampler.reset()

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