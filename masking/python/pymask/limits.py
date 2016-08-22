import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf
import pickle
import gzip
import pymultinest
import os, threading, subprocess
import matplotlib.pyplot as plt
import json
import pymask.oifits as oifits
import time

import pymask

'''---------------------------------------------------------
limits.py - a test script to make a contrast detection limit
finder in Python using pymask. This is a translation of Mike
Ireland's contrast detection limit code in binary_grid. 

For this first stab, data are assumed to be uncorrelated.
This will change in future versions.
---------------------------------------------------------'''

ddir = '/home/ben/Research/Desert/'

oifitsfile = 'Binary1.oifits' #'ScoCen2_Targ1.oifits'

#------------------------
# first, load your data!
#------------------------

a = pymask.cpo(ddir+oifitsfile)
t3err = a.t3err
ndata = a.ndata
u,v,wavel = a.u,a.v,a.wavel

#----------------------------
# Create your parameter space
#----------------------------


sepmin,sepmax = 20.0,200.0 # grid dimensions - angle is automatically
#from 0 to 360
initcon = 10.0 # pick an initial contrast ratio
nseps = 100 # number of points per axis in the grid
nangles = 100

seps = np.linspace(sepmin,sepmax,nseps)
angles = np.linspace(0,360,nangles,nangles)

params = np.zeros((3,nseps,nangles))

for i in range(nseps):
	for j in range(nangles):
		params[0,i,j] = seps[i]
		params[1,i,j] = angles[j]
		params[2,i,j] = initcon

#------------------------
# Monte Carlo time!
#------------------------


nsims = 10000

ncps = len(a.t3data)

tic = time.time()

#create random data
mcdata =np.random.randn(ncps,nsims)
for i in range(nsims):
	mcdata[:,i] *= t3err

maxrats = np.zeros((2*nseps,nsims)) # maximum contrast ratios for each separation
maxgtrats = np.zeros((2*nseps-1,nsims)) # maximum contrast ratios for separations outside of seps

#create a grid of models
print 'Initialising model grid'
modcp = np.zeros((ndata,nangles*nseps))
for i in range(nseps):
	for j in range(nangles):
		modcp[:, i+long(nseps)*j] = pymask.cp_model(params[:, i, j],a.u,a.v,a.wavel)

errormatrix = np.diag(1/t3err)

#create the likelihood norms
print 'Creating likelihood norms'
norm =  np.zeros((nseps*nangles)) 
for i in range(nseps*nangles):
	norm[i]=np.dot(np.dot(modcp[:,i],errormatrix),modcp[:,i])

norm = np.reshape(norm,(nseps,nangles))

# loop over the simulations!
print 'Looping over the simulations!'
for i in range(nsims):
	crats = 1./params[2, :, :]*np.reshape(np.dot(mcdata[:,i],errormatrix).dot(modcp), (nseps, nangles))/norm

crats = crats/((1.0 - 3.0*abs(crats))>0.4)

new_crats=pymask.rebin(crats,(2*nseps,2*nangles))

toc = time.time()

print 'Time elapsed =',toc-tic, 'seconds'

