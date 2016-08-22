#!/usr/bin/env python

import pysco
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pdb
import os

from numpy.random import rand, randn
from random import choice, shuffle

'''------------------------------------------------------------
Nested sampling code for kernel phase analysis.
------------------------------------------------------------'''

#ddir = './data/GO10879/' # data directory
ddir = "./data/GO10143/"

num = 10  # index of file in database

band = 2  # H = 1, J = 2

pupil = '14b'

#-----------------------------------------------------------------
# Initialise simulation parameters
#--------------------------------------------------------------'''

npoints = 300 # active points
nsteps = 6000 # iterations
njumps = 20 # random walk steps

sepmin,sepmax,anglemin,anglemax,cmin,cmax = 700.,750.,0.,360.,1.00001,2.0

# -------------------
# 1. load the dataset
# -------------------

if band==1:
    bandn = 'H'
elif band==2:
    bandn = 'J'
else: print 'Invalid band'

if num >= 10:
    num = str(num)
else:
    num = '0'+str(num)
    
if ddir == './data/GO10143/':
    fname = 'n8yj'+num+'0'+str(band)+'0_mos.fits'
elif ddir == './data/GO10879/':
    fname = 'n9nk'+num+'0'+str(band)+'0_mos.fits'
else: print 'Data directory failure'

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
a.name = "HST - " + str(num)+bandn #  # add a label to the template
a.extract_from_fits_frame(ddir+fname)
wavelength = a.info['filter']

a.load_cals('calibs'+pupil+'.pick')
if ddir == './data/GO10143/':
    if band == 1:
        a.calibrate(a.kmean10143H,a.krms10143H)
    elif band == 2:
        a.calibrate(a.kmean10143J,a.krms10143J)
    else: print 'Band failure'
elif ddir == './data/GO10879/':
    if band == 1:
        a.calibrate(a.kmean10879H,a.krms10879H)
    elif band == 2:
        a.calibrate(a.kmean10879J,a.krms10879J)
    else: print 'Band failure'
else: print 'Data directory failure'

#a.kp_signal = a.kp_error*np.random.randn(a.nkphi)
#shuffle(a.kp_signal)

logpre = a.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(a.kp_error))

def loglikelihood(sep,angle,contrast):
    '''Define the chi-squared function - omits logpre! '''
    
    test = phase_binary(a.uv[:,0],a.uv[:,1], wavelength,
                        [sep,angle+90.0-a.info['orient'],contrast])
    modl_ker = np.dot(a.KerPhi, a.RED*test)

    chisquared = np.sum(((modl_ker - a.kp_signal)/(a.kp_error))**2)

    #logpre = a.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(a.kp_error))
    #logpre = 0 
    #like = logpre-chisquared/2
    like = -chisquared/2
    
    return like


'''--------------------------------------------------------
Generate 100 active points from a uniform prior in angle,
and from a Jeffreys prior in separation and contrast
(IMPLEMENT!)
--------------------------------------------------------'''

seps = np.exp(rand(npoints)*np.log(sepmax/sepmin))*sepmin
angles = (anglemax-anglemin)*rand(npoints)+anglemin
cons = np.exp(rand(npoints)*np.log(cmax/cmin))*cmin
#cons = (cmax-cmin)*rand(npoints)*cmin

'''---------------------------------------------------------
Loop over 1000 iterations of the nested sampling algorithm.
---------------------------------------------------------'''

L = np.zeros(npoints)

rejected = np.zeros((nsteps,3))
rejL = np.zeros(nsteps)

lognull = -np.sum(((a.kp_signal)/a.kp_error)**2)/2 +logpre

progress = 99 # counter for displays

for k in range(0,npoints): # calculate likelihood

    L[k] = loglikelihood(seps[k],angles[k],cons[k])

for j in range(0,nsteps): # loop for iterations

    print 'Commencing iteration',j

    # store rejected binary parameters and re-sort L

    jumps = (sepmax-sepmin)/200. # step rescaling parameters
    jumpa = (anglemax-anglemin)/200.
    jumpc = (cmax-cmin)/200.
        
    rejected[j,0] = seps[L==L.min()]
    rejected[j,1] = angles[L==L.min()]
    rejected[j,2] = cons[L==L.min()]
    rejL[j] = L[L==L.min()]
    
    args = np.argsort(L)
    seps, angles, cons, L = seps[args], angles[args], cons[args], L[args]

    acceptances = 0 # initialise for coming loop
    rejections = 0
    l = 0
    trysd,tryad,trycd = 0,0,0

    choose = int(np.ceil(rand()*(npoints-1))) # make a random choice

    starts, starta, startc = seps[choose],angles[choose],cons[choose]
    news,newa,newc = starts, starta, startc
    newl = loglikelihood(starts,starta,startc)
    
    while l<njumps: # random walk to generate a new active point

        trys = news + jumps*randn()
        trya = np.mod(newa + jumpa*randn(),anglemax)
        tryc = newc + jumpc*randn()

        like = loglikelihood(trys,trya,tryc)
        
        if like > L[0] and sepmin<trys<sepmax and cmin<tryc<cmax: #accept and move on
            acceptances += 1
            news, newa, newc, newl = trys, trya, tryc, like
        else: # dwell on it
            rejections += 1

        #adjust jumps - factors from John Skilling's example

        if acceptances > rejections: # jump further if you can afford to
            jumps *= np.exp(1/acceptances)
            jumpa *= np.exp(1/acceptances)
            jumpc *= np.exp(1/acceptances)
            
        else: # don't jump as far if you're making it worse too often
            jumps /= np.exp(1/float(rejections))
            jumpa /= np.exp(1/float(rejections))
            jumpc /= np.exp(1/float(rejections))
            
        l += 1

        if (l >= njumps and acceptances ==0): # start somewhere else
            l=0
            choose = int(np.ceil(rand()*(npoints-1)))
            news,newa,newc = seps[choose],angles[choose],cons[choose]
            
    print l, 'steps taken', acceptances/float(l),'acceptance rate'
    print 'Sep',news,'Angle',newa,'Contrast',newc

    seps[0], angles[0], cons[0], L[0] = news, newa, newc, newl

    if progress == 99: # make a plot every 100 iterations
        print 'Generating plots'

        plt.figure(0)
        plt.clf()
        plt.scatter(seps,angles)
        plt.axis([sepmin,sepmax,anglemin,anglemax])
        plt.xlabel('Separation (mas)')
        plt.ylabel('Position Angle (degrees)')
        plt.title('Active Points')
        plt.draw()
        
        #plt.show()

        plt.figure(1)
        plt.clf()
        plt.scatter(seps,cons)
        plt.axis([sepmin,sepmax,cmin,cmax])
        plt.xlabel('Separation (mas)')
        plt.ylabel('Contrast Ratio')
        plt.title('Active Points')
        plt.draw()
        
        plt.show()
        
        progress = 0
    else:
        progress += 1

'''-----------------------------------------------
We can now integrate these results - we choose the
rectangle rule wi = xi-x(i-1)
-----------------------------------------------'''

x =  np.exp(-np.arange(0,nsteps)/float(npoints))

seps = rejected[:,0]
angles = rejected[:,1]
cons = rejected[:,2]

xx = np.zeros(len(x)+2) #dummy padded array for reflecting bcs
xx[0] = 2-x[0]
xx[-1]= x[-1]
xx[1:-1] = x

w = 0.5*np.abs(np.array([xx[j-1]-xx[j+1] for j in range(1,len(xx)-1)]))
logw = np.log(w)

E = np.sum(w*np.exp(rejL-rejL.max()))
logE = np.log(E) + rejL.max()

logpdf = rejL+logw-logE

logE += logpre # get the prefactor back in

pdf = np.exp(logpdf)
        
H = np.sum(pdf*(rejL-logE-logpre))

dlogE = np.sqrt(H/nsteps)

#dlnE = np.sqrt(H/nsteps)

sepmean = np.sum(pdf*seps)
dsep = np.sqrt(np.sum(pdf * (seps-sepmean)**2))

anglemean = np.sum(pdf*angles)
dangle = np.sqrt(np.sum(pdf*(angles-anglemean)**2))

cmean = np.sum(pdf*cons)
dc = np.sqrt(np.sum(pdf*(cons-cmean)**2))

'''----------------------------------------------
Return results and display correlation diagram
----------------------------------------------'''

print 'Log E =',logE+logpre,'dlogE =',dlogE, 'log null =',lognull,'significance =',logE+logpre-lognull

print 'Mean sep', sepmean, 'mean angle',anglemean,'mean contrast',cmean

print 'dsep',dsep,'dangle',dangle,'dc',dc

params = [sepmean,anglemean,cmean]
params2 = [seps[rejL==rejL.max()],angles[rejL==rejL.max()],cons[rejL==rejL.max()]]
b = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
b.binary_model(params, a.info)

#a.correlation_plot(b)

odds = np.exp(logE-lognull)
prob = 1./(1.+1./(odds))

print 'Odds over null hypothesis',odds
print 'Probability',prob

plt.clf()
mm = np.round(np.max(np.abs(a.kp_signal)), -1)
plt.errorbar(b.kp_signal,a.kp_signal,yerr=a.kp_error,fmt='b.')
plt.plot([-mm,mm],[-mm,mm], 'g')
plt.axis('tight',fontsize='large')
plt.xlabel('Model Kernel Phases',fontsize='large')
plt.ylabel('Kernel Phase Signal', fontsize='large')
plt.title('Kernel Phase Correlation Diagram, target '+str(num)+', '+str(bandn)+' band',
          fontsize='large')
plt.draw()
plt.show()

