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
from multiprocessing import Pool
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

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
	cps = np.array(np.sum(phases,axis=1))
	return cps

# =========================================================================
# =========================================================================

def hammer(cpo,ivar=[52., 192., 1.53],ndim=3,nwalcps=50,plot=False):

	'''Default implementation of emcee, the MCMC Hammer, for closure phase
	fitting. Requires a closure phase object cpo, and is best called with 
	ivar chosen to be near the peak - it can fail to converge otherwise.'''

	ivar = np.array(ivar)  # initial parameters for model-fit

	p0 = [ivar + 0.1*ivar*np.random.rand(ndim) for i in range(nwalcps)] # initialise walcps in a ball

	print 'Running emcee now!'

	t0 = time.time()

	sampler = emcee.EnsembleSampler(nwalcps, ndim, cp_loglikelihood, args=[cpo.u,cpo.v,cpo.wavel,cpo.t3data,cpo.t3err])
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

def nest(cpo,paramlimits=[20.,250.,0.,360.,1.0001,10],ndim=3,resume=False,eff=0.3,multi=True):

	'''Default implementation of a MultiNest fitting routine for closure 
	phase data. Requires a closure phase cpo object, parameter limits and 
	sensible keyword arguments for the multinest parameters. 

	This function does very naughty things creating functions inside this 
	function because PyMultiNest is very picky about how you pass it
	data.

	Optional parameter eff tunes sampling efficiency, and multi toggles multimodal 
	nested sampling on and off. Turning off multimodal sampling results in a speed 
	boost of ~ 20-30%. 

	'''

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

	tic = time.time() # start timing

	#---------------------------------
	# now run MultiNest!
	#---------------------------------

	pymultinest.run(myloglike, myprior, n_params, wrapped_params=[1], resume=resume, verbose=True, sampling_efficiency=eff, multimodal=multi, n_iter_before_update=1000)

	# let's analyse the results
	a = pymultinest.Analyzer(n_params = n_params)
	s = a.get_stats()

	toc = time.time()

	if toc-tic < 60.:
		print 'Time elapsed =',toc-tic,'s'
	else: 
		print 'Time elapsed =',(toc-tic)/60.,'mins'

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

# =========================================================================
# =========================================================================


def detec_sim_loopfit(everything):
	'''Function for multiprocessing in detec_limits. Takes a 
	single separation and full angle, contrast lists.'''
	chi2_diff = np.zeros((everything['nth'],everything['ncon'],everything['nsim']))
	for j,th in enumerate(everything['ths']):
		for k,con in enumerate(everything['cons']):
			bin_cp = cp_model([everything['sep'],th,con],everything['u'], everything['v'], everything['wavel'])
			
			# binary cp model
			# ----------------------
			
			rnd_cp = bin_cp[:,np.newaxis] + everything['error'][:,np.newaxis]*everything['rands']
			chi2_sngl = np.sum((((rnd_cp)/ everything['error'][:,np.newaxis])**2),axis=0) 
			chi2_binr = np.sum((((rnd_cp-bin_cp[:,np.newaxis]) / everything['error'][:,np.newaxis])**2),axis=0)
			chi2_diff[j,k] = chi2_binr-chi2_sngl# note not i,j,k - i is for seps
	if everything['ix'] % 50 ==0:
		print 'Done',everything['ix']
	return chi2_diff

# =========================================================================
# =========================================================================


def detec_limits(cpo,nsim=100,nsep=32,nth=20,ncon=32,smin='Default',smax='Default',
	cmin=1.0001,cmax=500.,addederror=0,threads=0,save=False):

	'''uses a Monte Carlo simulation to establish contrast-separation 
	detection limits given an array of standard deviations per closure phase.

	Because different separation-contrast grid points are entirely
	separate, this task is embarrassingly parallel. If you want to 
	speed up the calculation, use multiprocessing with a threads 
	argument equal to the number of available cores.

	Make nseps a multiple of threads! This uses the cores most efficiently.

	Hyperthreading (2x processes per core) in my experience gets a ~20%
	improvement in speed.

	Written by F. Martinache and B. Pope.'''

	#------------------------
	# first, load your data!
	#------------------------

	error = cpo.t3err + addederror

	u,v = cpo.u,cpo.v

	wavel = cpo.wavel

	ndata = cpo.ndata

	w = np.array(np.sqrt(u**2 + v**2))/wavel

	if smin == 'Default':
		smin = rad2mas(1./4/np.max(w))

	if smax == 'Default':
		smax = rad2mas(1./np.min(w))

	#------------------------
	# initialise Monte Carlo
	#------------------------

	seps = smin + (smax-smin) * np.linspace(0,1,nsep)
	ths  = 0.0 + 360.0  * np.linspace(0,1,nth)
	cons = cmin  + (cmax-cmin)  * np.linspace(0,1,ncon)

	rands = np.random.randn(ndata,nsim)

	#------------------------
	# Run Monte Carlo
	#------------------------

	tic = time.time() # start the clock

	if threads ==0:
		chi2_diff = np.zeros((nsep,nth,ncon,nsim))
		for i,sep in enumerate(seps):
			print("iteration # %3d: sep=%.2f" % (i, sep))
			chi2_diff[i,:,:,:]= detec_sim_loopfit(sep)
			toc = time.time()
			if i != 0:
				remaining =  (toc-tic)*(nsep-i)/float(i)
				if remaining > 60:
					print('Estimated time remaining: %.2f mins' % (remaining/60.))
				else: 
					print('Estimated time remaining: %.2f seconds' % (remaining))
	else:
		all_vars=[]
		for ix in range(nsep):
			everything={'sep':seps[ix],'cons':cons,'ths':ths, 'ix':ix,
				'nsep':nsep,'ncon':ncon,'nth':nth,'u':u,'v':v,'nsim':nsim,
				'rands':rands,'error':error,'wavel':wavel}
			all_vars.append(everything)
		pool = Pool(processes=threads)
		chi2_diff=pool.map(detec_sim_loopfit,all_vars)
		chi2_diff = np.array(chi2_diff)
	tf = time.time()
	if tf-tic > 60:
		print 'Total time elapsed:',(tf-tic)/60.,'mins'
	elif tf-tic <= 60:
		print 'Total time elapsed:',tf-tic,'seconds'

	ndetec = np.zeros((ncon, nsep))

	nc, ns = int(ncon), int(nsep)

	for k in range(nc):
		for i in range(ns):
			toto = (chi2_diff)[i,:,k,:]
			ndetec[k,i] = (toto < 0.0).sum()

	nbtot = nsim * nth

	ndetec /= float(nbtot)

	# ---------------------------------------------------------------
	#                        contour plot!
	# ---------------------------------------------------------------
	levels = [0.5,0.9, 0.99, 0.999]
	mycols = ('k', 'k', 'k', 'k')

	plt.figure(0)
	contours = plt.contour(ndetec, levels, colors=mycols, linewidth=2, 
                 extent=[smin, smax, cmin, cmax])
	plt.clabel(contours)
	plt.contourf(seps,cons,ndetec,levels,cmap=plt.cm.bone)
	plt.colorbar()
	plt.xlabel('Separation (mas)')
	plt.ylabel('Contrast Ratio')
	plt.title('Contrast Detection Limits')
	plt.draw()
	plt.show()

	data = {'levels': levels,
			'ndetec': ndetec,
			'seps'  : seps,
			'angles': ths,
			'cons'  : cons,
			'name'  : cpo.name}

	if save == True:
		file = 'limit_lowc'+cpo.name+'.pick'
		print file

		myf = open(file,'w')
		pickle.dump(data,myf)
		myf.close()

	return data

# =========================================================================
# =========================================================================

def binary_fit(cpo, p0):
    '''Performs a best binary fit search for the dataset.
    -------------------------------------------------------------
    p0 is the initial guess for the parameters 3 parameter vector
    typical example would be : [100.0, 0.0, 5.0].
    returns the full solution of the least square fit:
    - soluce[0] : best-fit parameters
    - soluce[1] : covariance matrix
    ------------------------------------------------------------- '''
    
    if np.all(cpo.t3err == 0.0):
        print("Closure phase object instance is not calibrated.\n")
        soluce = leastsq(cpo.bin_fit_residuals, p0, args=(cpo),
                     full_output=1)
    else:
    	def lmcpmodel(index,params1,params2,params3):
    		params = [params1,params2,params3]
    		model = cp_model(params,cpo.u,cpo.v,cpo.wavel)
    		return model[index]
    	soluce = curve_fit(lmcpmodel,range(0,cpo.ndata),cpo.t3data,p0,sigma=cpo.t3err)
    cpo.covar = soluce[1]
    soluce[0][1] = np.mod(soluce[0][1],360.) # to get consistent position angle measurements
    return soluce

# =========================================================================
# =========================================================================

def bin_fit_residuals(params, cpo):
	'''Function for binary_fit without errorbars'''
	test = cp_model(params,cpo.u,cpo.v,cpo.wavel)
	err = (cpo.t3data - test)
	return err
	
# =========================================================================
# =========================================================================
	
def brute_force_chi2_grid(everything):
	'''Function for multiprocessing, does 3d chi2 fit, followed by
	   Levenberg Marquadt, then returns best sep, PA, contrast ratio'''
	this_cpo=everything['sim_cpo']
	data_cp=this_cpo.t3data
	chi2=np.zeros((everything['nsep'],everything['nth'],everything['ncon']))
	for i,sep in enumerate(everything['seps']):
		for j,th in enumerate(everything['ths']):
			for k,con in enumerate(everything['cons']):
				mod_cps = cp_model([sep,th,con],everything['u'],everything['v'],everything['wavel'])
				chi2[i,j,k]=np.sum(((data_cp-mod_cps)/everything['error'])**2)
	b_params_ix=np.where(chi2==np.amin(chi2))
	b_params=[everything['seps'][b_params_ix[0][0]],everything['ths'][b_params_ix[1][0]],everything['cons'][b_params_ix[2][0]]]
	b_params=np.array(b_params)
	
	#now do L-M fitting. Sometimes it can't find the best position,
	#so take the coarse estimate instead
	try:
		[best_params,cov]=binary_fit(this_cpo,b_params)
	except:
		best_params=b_params
		print "couldn't find best params!"
	if everything['ix'] % 50 ==0:
		print 'Done',everything['ix']
	cov,chi2,b_params,this_cpo,data_cp=None,None,None,None,None
	return best_params

# =========================================================================
# =========================================================================


def brute_force_detec_limits(cpo,nsim=100,nsep=32,nth=20,ncon=32,smin='Default',smax='Default',
	cmin=10.,cmax=500.,addederror=0,threads=0,save=False):

	'''uses a Monte Carlo simulation to establish contrast-separation 
	detection limits given an array of standard deviations per closure phase.

	Because different separation-contrast grid points are entirely
	separate, this task is embarrassingly parallel. If you want to 
	speed up the calculation, use multiprocessing with a threads 
	argument equal to the number of available cores.

	Make nseps a multiple of threads! This uses the cores most efficiently.

	Hyperthreading (2x processes per core) in my experience gets a ~20%
	improvement in speed.

	Written by F. Martinache and B. Pope.
	
	This version was modified by ACC to use a brute force
	chi2 grid instead of making any assumptions or approximations (and is
	obviously slower).'''

	#------------------------
	# first, load your data!
	#------------------------

	error = cpo.t3err + addederror

	u,v = cpo.u,cpo.v

	wavel = cpo.wavel

	ndata = cpo.ndata

	w = np.array(np.sqrt(u**2 + v**2))/wavel

	if smin == 'Default':
		smin = rad2mas(1./4/np.max(w))

	if smax == 'Default':
		smax = rad2mas(1./np.min(w))

	#------------------------
	# initialise Monte Carlo
	#------------------------

	seps = smin + (smax-smin) * np.linspace(0,1,nsep)
	ths  = 0.0 + 360.0  * np.linspace(0,1,nth)
	cons = cmin  + (cmax-cmin)  * np.linspace(0,1,ncon)

	rands = np.random.randn(ndata,nsim)
	sim_cps=rands*error[:,np.newaxis]
	sim_cpos=[]
	for ix in range(nsim):
		this_cpo=copy.deepcopy(cpo)
		this_cpo.t3data=sim_cps[:,ix]
		sim_cpos.append(this_cpo)
	
	#------------------------
	# Run Monte Carlo
	#------------------------

	tic = time.time() # start the clock
	best_params=[]
	if threads ==0:
		toc=time.time()
		best_params=np.zeros((nsim,3))
		for ix in range(nsim):
			best_params[ix,:]=brute_force_chi2_grid(ix)
			if (ix % 50) ==0:
				tc=time.time()
				print 'Done',ix,'. Time taken:',(tc-toc),'seconds'
				toc=tc	
	else:
		all_vars=[]
		for ix in range(nsim):
			everything={'seps':seps,'cons':cons,'ths':ths, 'ix':ix,
				'nsep':nsep,'ncon':ncon,'nth':nth,'u':u,'v':v,
				'sim_cpo':sim_cpos[ix],'error':error,'wavel':wavel}
			all_vars.append(everything)
		pool = Pool(processes=threads)
		best_params=pool.map(brute_force_chi2_grid,all_vars)
	tf = time.time()
	if tf-tic > 60:
		print 'Total time elapsed:',(tf-tic)/60.,'mins'
	elif tf-tic <= 60:
		print 'Total time elapsed:',tf-tic,'seconds'

	ndetec = np.zeros((ncon, nsep))
	nc, ns = int(ncon), int(nsep)
	
	#collect them	
	for ix in range(nsim):
		sep=best_params[ix][0]
		con=best_params[ix][2]
		sep_ix=np.where(abs(seps-sep)==np.amin(abs(seps-sep)))
		con_ix=np.where(abs(cons-con)==np.amin(abs(cons-con)))
		ndetec[sep_ix,con_ix]+=1
	#Take the cumulative sum over contrast ratio at each sep
	cumsum_detec=ndetec.cumsum(axis=1)
	#turn into %
	maxdetec=np.amax(cumsum_detec,axis=1)
	ndetec=0*cumsum_detec
	for ix in range(nsep):
		if maxdetec[ix]==0:
			print 'No sims for sep '+str(seps[ix])+'mas'
		else:
			print str(maxdetec[ix])+" in "+str(seps[ix])+" bin."
			ndetec[ix,:]=cumsum_detec[ix,:]/maxdetec[ix]

	ndetec=1-ndetec
	#Axes were wrong way around (I blame IDL)
	ndetec=np.transpose(ndetec)
	
	# ---------------------------------------------------------------
	#                        contour plot!
	# ---------------------------------------------------------------
	levels = [0.,0.9, 0.99, 0.999]
	mycols = ('k', 'k', 'k', 'k')

	plt.figure(0)
	contours = plt.contour(ndetec, levels, colors=mycols, linewidth=2, 
		 extent=[smin, smax, cmin, cmax])
	plt.clabel(contours)
	plt.contourf(seps,cons,ndetec,levels,cmap=plt.cm.bone)
	plt.colorbar()
	plt.xlabel('Separation (mas)')
	plt.ylabel('Contrast Ratio')
	plt.title('Contrast Detection Limits')
	plt.draw()
	plt.show()

	data = {'levels': levels,
			'ndetec': ndetec,
			'seps'  : seps,
			'angles': ths,
			'cons'  : cons,
			'name'  : cpo.name,
			'cumsum': cumsum_detec,
			'best_params': best_params}

	if save == True:
		file = 'limit_lowc'+cpo.name+'.pick'
		print file

		myf = open(file,'w')
		pickle.dump(data,myf)
		myf.close()
		
	return data