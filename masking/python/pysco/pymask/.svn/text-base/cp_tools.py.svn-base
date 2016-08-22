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
import os, threading, subprocess
import matplotlib.pyplot as plt
import json
import oifits
import time
import emcee

'''------------------------------------------------------------------------
cp_tools.py - a collection of functions useful for closure phase analysis
in Python. This includes mas2rad, rad2mas and phase_binary from pysco; 
it must share a directory with oifits.py, and depends on PyMultiNest,
MultiNest (Fortran) and emcee (the Python MCMC Hammer implementation).
------------------------------------------------------------------------'''

def mas2rad(x):
    ''' Convenient little function to convert milliarcsec to radians '''
    return x*np.pi/(180*3600*1000)

# =========================================================================
# =========================================================================

def rad2mas(x):
    ''' Convenient little function to convert radians to milliarcseconds '''
    return x/np.pi*(180*3600*1000)
# =========================================================================
# =========================================================================

def phase_binary(u, v, wavel, p):
    ''' Calculate the phases observed by an array on a binary star
    ----------------------------------------------------------------
    p: 3-component vector (+2 optional), the binary "parameters":
    - p[0] = sep (mas)
    - p[1] = PA (deg) E of N.
    - p[2] = contrast ratio (primary/secondary)
    
    optional:
    - p[3] = angular size of primary (mas)
    - p[4] = angular size of secondary (mas)

    - u,v: baseline coordinates (meters)
    - wavel: wavelength (meters)
    ---------------------------------------------------------------- '''

    p = np.array(p)
    # relative locations
    th = (p[1] + 90.0) * np.pi / 180.0
    ddec =  mas2rad(p[0] * np.sin(th))
    dra  = -mas2rad(p[0] * np.cos(th))

    # baselines into number of wavelength
    x = np.sqrt(u*u+v*v)/wavel

    # decompose into two "luminosity"
    l2 = 1. / (p[2] + 1)
    l1 = 1 - l2
    
    # phase-factor
    phi = np.zeros(u.size, dtype=complex)
    phi.real = np.cos(-2*np.pi*(u*dra + v*ddec)/wavel)
    phi.imag = np.sin(-2*np.pi*(u*dra + v*ddec)/wavel)

    # optional effect of resolved individual sources
    if p.size == 5:
        th1, th2 = mas2rad(p[3]), mas2rad(p[4])
        v1 = 2*j1(np.pi*th1*x)/(np.pi*th1*x)
        v2 = 2*j1(np.pi*th2*x)/(np.pi*th2*x)
    else:
        v1 = np.ones(u.size)
        v2 = np.ones(u.size)

    cvis = l1 * v1 + l2 * v2 * phi
    phase = np.angle(cvis, deg=True)
    return np.mod(phase + 10980., 360.) - 180.0

# =========================================================================
# =========================================================================

def cp_loglikelihood(params,u,v,wavel,t3data,t3err):
    '''Calculate loglikelihood for closure phase data.
    Used both in the MultiNest and MCMC Hammer implementations.'''
    cps = cp_model(params,u,v,wavel)
    chi2 = np.sum(((t3data-cps)/t3err)**2)
    loglike = -chi2/2
    return loglike

# =========================================================================
# =========================================================================

def cp_model(params,u,v,wavel):
    '''Function to model closure phases. Takes a parameter list, u,v triangles and a single wavelength.'''
    ndata = u.shape[0]
    phases = phase_binary(u.ravel(),v.ravel(),wavel,params)
    phases = np.reshape(phases,(ndata,3))
    cps = np.sum(phases,axis=1)
    return cps

# =========================================================================
# =========================================================================

def hammer(cpo,ivar=[52., 192., 1.53],ndim=3,nwalkers=50,plot=False):

    '''Default implementation of emcee, the MCMC Hammer, for closure phase
    fitting. Requires a closure phase object cpo, and is best called with 
    ivar chosen to be near the peak - it can fail to converge otherwise.'''

    ivar = np.array(ivar)  # initial parameters for model-fit

    p0 = [ivar + 0.1*ivar*np.random.rand(ndim) for i in range(nwalkers)] # initialise walkers in a ball

    print 'Running emcee now!'

    t0 = time.time()

    sampler = emcee.EnsembleSampler(nwalkers, ndim, cp_loglikelihood, args=[cpo.u,cpo.v,cpo.wavel,cpo.t3data,cpo.t3err])
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

    if plot==True:

        plt.clf()

        paramnames = ['Separation','Position Angle','Contrast']
        paramdims = ['(mas)', '(deg)','Ratio']

        for i in range(ndim):
            plt.figure(i)
            plt.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
            plt.title(paramnames[i])
            plt.ylabel('Counts')
            plt.xlabel(paramnames[i]+paramdims[i])

        plt.show()

    # =========================================================================
# =========================================================================

def nest(cpo,paramlimits=[20.,250.,0.,360.,1.0001,10],ndim=3,resume=False,eff=0.3):

    '''Default implementation of a MultiNest fitting routine for closure 
    phase data. Requires a closure phase cpo object, parameter limits and 
    sensible keyword arguments for the multinest parameters. 

    This function does very naughty things creating functions inside this 
    function because PyMultiNest is very picky about how you pass it
    data.'''

    def myprior(cube, ndim, n_params,paramlimits=paramlimits,kpo=0):
        cube[0] *= (paramlimits[1] - paramlimits[0])+paramlimits[0]
        cube[1] *= (paramlimits[3] - paramlimits[2])+paramlimits[2]
        cube[2] *= (paramlimits[5] - paramlimits[4])+paramlimits[4]

    def myloglike(cube, ndim, n_params):
        loglike = cp_loglikelihood(cube[0:3],cpo.u,cpo.v,cpo.wavel,cpo.t3data,cpo.t3err)
        return loglike

    parameters = ['Separation','Position Angle','Contrast']
    n_params = len(parameters)
    ndim = 3

    pymultinest.run(myloglike, myprior,n_params, wrapped_params=[1], resume = resume, verbose = True, sampling_efficiency = eff)

    # lets analyse the results
    a = pymultinest.Analyzer(n_params = n_params)
    s = a.get_stats()

    # json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
    print
    print "-" * 30, 'ANALYSIS', "-" * 30
    print "Global Evidence:\n\t%.15e +- %.15e" % ( s['global evidence'], s['global evidence error'] )
    params = s['marginals']

    bestsep = params[0]['median']
    seperr = params[0]['sigma']

    bestth = params[1]['median']
    therr = params[1]['sigma']

    bestcon = params[2]['median']
    conerr = params[2]['sigma']
    
    print ''

    print 'Separation',bestsep,'pm',seperr
    print 'Position angle',bestth,'pm',therr
    print 'Contrast ratio',bestcon,'pm',conerr

    return s