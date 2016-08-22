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
import pymultinest
import math
import os, threading, subprocess
import matplotlib.pyplot as plt
if not os.path.exists("chains"): os.mkdir("chains")
import json
import oifits
import time
import pysco

'''---------------------------------------------------
Script to try out PyMultiNest 
---------------------------------------------------'''
# our probability functions
# Taken from the eggbox problem.

def cp_model(params,u,v,wavel):
    '''Function to model closure phases. Takes a parameter list, u,v triangles and a single wavelength.'''
    ndata = u.shape[0]
    phases = pysco.phase_binary(u.ravel(),v.ravel(),wavel,params)
    phases = np.reshape(phases,(ndata,3))
    cps = np.sum(phases,axis=1)
    return cps

def myprior(cube, ndim, n_params,paramlimits=[20.,250.,0.,360.,1.0001,10]):
	cube[0] *= (paramlimits[1] - paramlimits[0])+paramlimits[0]
	cube[1] *= (paramlimits[3] - paramlimits[2])+paramlimits[2]
	cube[2] *= (paramlimits[5] - paramlimits[4])+paramlimits[4]

'''----------------------------------------------------
First, read in the data
----------------------------------------------------'''

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

def loglikelihood(params,u,v,wavel,t3data,t3err):
    cps = cp_model(params,u,v,wavel)
    chi2 = np.sum(((t3data-cps)/t3err)**2)
    loglike = -chi2/2
    return loglike

def myloglike(cube, ndim, n_params):
	loglike = loglikelihood(cube[0:3],u,v,wavel,t3data,t3err)
	return loglike


'''----------------------------------------------------
Now we set up the multinest routine 
----------------------------------------------------'''
# number of dimensions our problem has
parameters = ['Separation','Position Angle','Contrast']
n_params = len(parameters)
ndim = 3

pymultinest.run(myloglike, myprior,n_params, wrapped_params=[1], resume = False, verbose = True, sampling_efficiency = 0.3)

# lets analyse the results
a = pymultinest.Analyzer(n_params = n_params)
s = a.get_stats()

# json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
print
print "-" * 30, 'ANALYSIS', "-" * 30
print "Global Evidence:\n\t%.15e +- %.15e" % ( s['global evidence'], s['global evidence error'] )

plt.clf()


