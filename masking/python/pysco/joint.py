import numpy as np
import matplotlib.pyplot as plt
import pyker
import pickle
from pyker_tools import *

from numpy.random import rand, randn
from random import choice, shuffle

'''------------------------------------------------------------
Nested sampling code for kernel phase analysis.
------------------------------------------------------------'''

ddir = './data/GO10879/' # data directory
#ddir = "./data/GO10143/"

num = 21 # index of file in database

pupil = '12b'
name = '2028+0052'

params = [45.8,107.7,3.1,1.52]

addederror = [0,0]

#-----------------------------------------------------------------
# Initialise simulation parameters
#--------------------------------------------------------------'''

npoints = 250 # active points
nsteps = 7500 # iterations
njumps = 20 # random walk steps

sepmin,sepmax,anglemin,anglemax,cmin,cmax = 30.,160.,0.,360.,1.001,3.5


# -------------------
# 1. load the dataset
# -------------------

if num >= 10:
    num = str(num)
else:
    num = '0'+str(num)
    
if ddir == './data/GO10143/':
    fnameH = 'n8yj'+num+'0'+str(1)+'0_mos.fits'
    fnameJ = 'n8yj'+num+'0'+str(2)+'0_mos.fits'
elif ddir == './data/GO10879/':
    fnameH = 'n9nk'+num+'0'+str(1)+'0_mos.fits'
    fnameJ = 'n9nk'+num+'0'+str(2)+'0_mos.fits'
elif ddir == './data/tdwarf/':
    fnameH = 'n8nh'+num+'0'+str(2)+'0_mos.fits'
    fnameJ = 'n8nh'+num+'0'+str(1)+'0_mos.fits'
else: print 'Data directory failure'

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
a.name = "HST - " + str(num)+'H' #  # add a label to the template
a.extract_from_fits_frame(ddir+fnameH)
wavelengthh = a.info['filter']

b = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
b.name = "HST - " + str(num)+'J' #  # add a label to the template
b.extract_from_fits_frame(ddir+fnameJ)
wavelengthj = b.info['filter']


if ddir == './data/GO10143/':
    a.load_cals('calibs'+pupil+'.pick')
    a.calibrate(a.kmean10143H,a.krms10143H)
    b.calibrate(a.kmean10143J,a.krms10143J)
elif ddir == './data/GO10879/':
    a.load_cals('calibs'+pupil+'.pick')
    a.calibrate(a.kmean10879H,a.krms10879H)
    b.calibrate(a.kmean10879J,a.krms10879J)
    
else: print 'Data directory failure'
#a.kp_signal = a.kp_error*np.random.randn(a.nkphi)
#shuffle(a.kp_signal)

##logpreh = a.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(a.kp_error))
##logprej = b.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(b.kp_error))
##
##def loglikelihood(sep,angle,contrasth,contrastj):
##    '''Define the chi-squared function - omits logpre! '''
##    
##    testh = phase_binary(a.uv[:,0],a.uv[:,1], wavelengthh,
##                        [sep,angle+90.0-b.info['orient'],contrasth])
##    modl_kerh = np.dot(a.KerPhi, a.RED*testh)
##
##    chisquaredh = np.sum(((modl_kerh - a.kp_signal)/(a.kp_error+addederror[0]))**2)
##
##    testj = phase_binary(b.uv[:,0],b.uv[:,1], wavelengthj,
##                        [sep,angle+90.0-b.info['orient'],contrastj])
##    modl_kerj = np.dot(b.KerPhi, b.RED*testj)
##
##    chisquaredj = np.sum(((modl_kerj - b.kp_signal)/(b.kp_error+addederror[1]))**2)
##
##    #logpre = a.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(a.kp_error))
##    #logpre = 0 
##    #like = logpre-chisquared/2
##
##    chisquared = chisquaredh + chisquaredj
##    
##    like = -chisquared/2
##    
##    return like
##
##
##'''--------------------------------------------------------
##Generate 100 active points from a uniform prior in angle,
##and from a Jeffreys prior in separation and contrast
##(IMPLEMENT!)
##--------------------------------------------------------'''
##
##seps = np.exp(rand(npoints)*np.log(sepmax/sepmin))*sepmin
##angles = (anglemax-anglemin)*rand(npoints)+anglemin
##consh = np.exp(rand(npoints)*np.log(cmax/cmin))*cmin
##consj = np.exp(rand(npoints)*np.log(cmax/cmin))*cmin
###cons = (cmax-cmin)*rand(npoints)*cmin
##
##'''---------------------------------------------------------
##Loop over 1000 iterations of the nested sampling algorithm.
##---------------------------------------------------------'''
##
##L = np.zeros(npoints)
##
##rejected = np.zeros((nsteps,4))
##rejL = np.zeros(nsteps)
##
##nullchi2 = np.sum(((a.kp_signal/(a.kp_error+addederror[0]))**2) + np.sum(((b.kp_signal)/(b.kp_error+addederror[1])))**2)
##lognull = -nullchi2/2 +logpreh + logprej
##
##progress = 99 # counter for displays
##
##for k in range(0,npoints): # calculate likelihood
##
##    L[k] = loglikelihood(seps[k],angles[k],consh[k],consj[k])
##
##for j in range(0,nsteps): # loop for iterations
##
##    print 'Commencing iteration',j
##
##    # store rejected binary parameters and re-sort L
##
##    jumps = (sepmax-sepmin)/200. # step rescaling parameters
##    jumpa = (anglemax-anglemin)/200.
##    jumpch = (cmax-cmin)/200.
##    jumpcj = (cmax-cmin)/200.
##        
##    rejected[j,0] = seps[L==L.min()]
##    rejected[j,1] = angles[L==L.min()]
##    rejected[j,2] = consh[L==L.min()]
##    rejected[j,3] = consj[L==L.min()]
##    rejL[j] = L[L==L.min()]
##    
##    args = np.argsort(L)
##    seps, angles, consh,consj, L = seps[args], angles[args], consh[args], consj[args],L[args]
##
##    acceptances = 0 # initialise for coming loop
##    rejections = 0
##    l = 0
##    trysd,tryad,trychd,trycjd = 0,0,0,0
##
##    choose = int(np.ceil(rand()*(npoints-1))) # make a random choice
##
##    starts, starta, startch,startcj = seps[choose],angles[choose],consh[choose],consj[choose]
##    news,newa,newch,newcj = starts, starta, startch,startcj
##    newl = loglikelihood(starts,starta,startch,startcj)
##    
##    while l<njumps: # random walk to generate a new active point
##
##        trys = news + jumps*randn()
##        trya = np.mod(newa + jumpa*randn(),anglemax)
##        trych = newch + jumpch*randn()
##        trycj = newcj + jumpcj*randn()
##
##        like = loglikelihood(trys,trya,trych,trycj)
##        
##        if like > L[0] and sepmin<trys<sepmax and cmin<trych<cmax and cmin<trycj<cmax and anglemin<trya<anglemax: #accept and move on
##            acceptances += 1
##            news, newa, newch, newcj, newl = trys, trya, trych, trycj, like
##        else: # dwell on it
##            rejections += 1
##
##        #adjust jumps - factors from John Skilling's example
##
##        if acceptances > rejections: # jump further if you can afford to
##            jumps *= np.exp(1/acceptances)
##            jumpa *= np.exp(1/acceptances)
##            jumpch *= np.exp(1/acceptances)
##            jumpch *= np.exp(1/acceptances)
##            
##        else: # don't jump as far if you're making it worse too often
##            jumps /= np.exp(1/float(rejections))
##            jumpa /= np.exp(1/float(rejections))
##            jumpch /= np.exp(1/float(rejections))
##            jumpcj /= np.exp(1/float(rejections))
##            
##        l += 1
##
##        if (l >= njumps and acceptances ==0): # start somewhere else
##            l=0
##            choose = int(np.ceil(rand()*(npoints-1)))
##            news,newa,newch,newcj = seps[choose],angles[choose],consh[choose],consj[choose]
##            
##    print l, 'steps taken', acceptances/float(l),'acceptance rate'
##    print 'Sep',news,'Angle',newa,'H band contrast',newch, 'J band contrast',newcj
##
##    seps[0], angles[0], consh[0],consj[0], L[0] = news, newa, newch, newcj, newl
##
##    if progress == 99: # make a plot every 100 iterations
##        print 'Generating plots'
##
##        plt.figure(0)
##        plt.clf()
##        plt.scatter(seps,angles)
##        plt.axis([sepmin,sepmax,anglemin,anglemax])
##        plt.xlabel('Separation (mas)')
##        plt.ylabel('Position Angle (degrees)')
##        plt.title('Active Points')
##        plt.draw()
##        
##        #plt.show()
##
##        plt.figure(1)
##        plt.clf()
##        plt.scatter(seps,consh)
##        plt.axis([sepmin,sepmax,cmin,cmax])
##        plt.xlabel('Separation (mas)')
##        plt.ylabel('Contrast Ratio (H Band)')
##        plt.title('Active Points')
##        plt.draw()
##
##        plt.figure(2)
##        plt.clf()
##        plt.scatter(consh,consj)
##        plt.axis([cmin,cmax,cmin,cmax])
##        plt.xlabel('Contrast Ratio (H Band)')
##        plt.ylabel('Contrast Ratio (J Band)')
##        plt.title('Active Points')
##        plt.draw()
##        
##        plt.show()
##        
##        progress = 0
##    else:
##        progress += 1
##
##'''-----------------------------------------------
##We can now integrate these results - we choose the
##rectangle rule wi = xi-x(i-1)
##-----------------------------------------------'''
##
##x =  np.exp(-np.arange(0,nsteps)/float(npoints))
##
##seps = rejected[:,0]
##angles = rejected[:,1]
##consh = rejected[:,2]
##consj = rejected[:,3]
##
##xx = np.zeros(len(x)+2) #dummy padded array for reflecting bcs
##xx[0] = 2-x[0]
##xx[-1]= x[-1]
##xx[1:-1] = x
##
##w = 0.5*np.abs(np.array([xx[j-1]-xx[j+1] for j in range(1,len(xx)-1)]))
##logw = np.log(w)
##
##E = np.sum(w*np.exp(rejL-rejL.max()))
##logE = np.log(E) + rejL.max()
##
##logpdf = rejL+logw-logE
##logpdf = logpdf-logpdf.max()
##
##pdf = np.exp(logpdf)
##pdf /= np.sum(pdf)
##
##logE += logpreh+logprej # get the prefactor back in
##
##H = np.sum(pdf*(rejL-logE-logpreh-logprej))
##
##dlogE = np.sqrt(H/nsteps)
##
###dlnE = np.sqrt(H/nsteps)
##
##sepmean = np.sum(pdf*seps)
##dsep = np.sqrt(np.sum(pdf * (seps-sepmean)**2))
##
##anglemean = np.sum(pdf*angles)
##dangle = np.sqrt(np.sum(pdf*(angles-anglemean)**2))
##
##chmean = np.sum(pdf*consh)
##dch = np.sqrt(np.sum(pdf*(consh-chmean)**2))
##
##cjmean = np.sum(pdf*consj)
##dcj = np.sqrt(np.sum(pdf*(consj-cjmean)**2))
##
##'''----------------------------------------------
##Return results and display correlation diagram
##----------------------------------------------'''
##
##print 'Log E =',logE+logpreh+logprej,'dlogE =',dlogE, 'log null =',lognull,'significance =',logE+logpreh+logprej-lognull
##
##print 'Mean sep', sepmean, 'mean angle',anglemean,'mean H contrast',chmean, 'mean J contrast',cjmean
##
##print 'dsep',dsep,'dangle',dangle,'dch',dch,'dcj',dcj
##
##params = [sepmean,anglemean,chmean,cjmean]
##uncertainties = [dsep,dangle,dch,dcj]
##
##paramsh = [sepmean,anglemean,chmean]
##paramsj = [sepmean,anglemean,cjmean]
##
paramsh= [params[0],params[1],params[2]]
paramsj= [params[0],params[1],params[3]]

#params2 = [seps[rejL==rejL.max()],angles[rejL==rejL.max()],cons[rejL==rejL.max()]]
c = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
c.binary_model(paramsh, a.info)

#a.correlation_plot(c)

d = pyker.KernelPhaseData("./kerphi/hst"+pupil+".pick")
d.binary_model(paramsj, b.info)

#b.correlation_plot(d)
##
##odds = np.exp(logE-lognull)
##prob = 1./(1.+1./(odds))
##
##print 'Odds over null hypothesis',odds
##print 'Probability',prob

plt.figure(3)
plt.clf()
mm = np.round(np.max(np.abs(a.kp_signal)), -1)
plt.errorbar(c.kp_signal,a.kp_signal,yerr=a.kp_error+addederror[0],fmt='b.')
plt.plot([-mm,mm],[-mm,mm], 'g')
plt.axis('tight',fontsize='large')
plt.xlabel('Model Kernel Phases',fontsize='large')
plt.ylabel('Kernel Phase Signal', fontsize='large')
plt.title('Kernel Phase Correlation Diagram, target '+name+' , H band',
          fontsize='large')
plt.draw()

plt.figure(4)
plt.clf()
mm = np.round(np.max(np.abs(b.kp_signal)), -1)
plt.errorbar(d.kp_signal,b.kp_signal,yerr=b.kp_error+addederror[1],fmt='b.')
plt.plot([-mm,mm],[-mm,mm], 'g')
plt.axis('tight',fontsize='large')
plt.xlabel('Model Kernel Phases',fontsize='large')
plt.ylabel('Kernel Phase Signal', fontsize='large')
plt.title('Kernel Phase Correlation Diagram, target '+name+' , J band',
          fontsize='large')
plt.draw()
plt.show()


