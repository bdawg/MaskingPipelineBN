import numpy as np
from scipy.optimize import leastsq
from core import *
import sys
from numpy.random import rand, randn
from random import choice, shuffle

# =========================================================================
# =========================================================================
def binary_model(params, kpi, hdr):
    ''' Creates a binary Kernel-phase model.
    
    ------------------------------------------------------------------ 
    uses a simple 5 parameter binary star model for the uv phases that
    should be observed with the provided geometry.
    
    Additional parameters are:
    - kpi, a kernel phase info structure
    - hdr, a header data information
    ------------------------------------------------------------------ '''
    
    params2 = np.copy(params)
    if 'Hale' in hdr['tel']: params2[1] += 220.0 + hdr['orient']
    if 'HST'  in hdr['tel']: params2[1] -= hdr['orient']
    else:                    params2[1] += 0.0

    filter = hdr['filter']
    
    testPhi = phase_binary(kpi.uv[:,0], kpi.uv[:,1], filter, params2)
    
    kpm = np.dot(kpi.KerPhi, kpi.RED * testPhi)
    return kpm

# =========================================================================
# =========================================================================
def binary_KPD_model(kpo, params):
    ''' Returns a 1D or 2D array of binary models.
    
    parameters are:
    - the kp structure (kpo) to be used as a template
    - the list of parameters of the binary model. '''

    models = np.zeros_like(kpo.kpd)

    nm = np.size(kpo.hdr) # number of kp realizations in kpo

    if nm == 1:
        models = binary_model(params, kpo.kpi, kpo.hdr)
    else:
        for i in xrange(nm):
            models[i] = binary_model(params, kpo.kpi, kpo.hdr[i])

    return models

# =========================================================================
# =========================================================================
def binary_KPD_fit_residuals(params, kpo):
    ''' Function to evaluate fit residuals, to be used in a leastsq
    fitting procedure. '''

    #test = binary_model(params, kpo.kpi, kpo.hdr)
    #test = binary_model(params, kpo.kpi, kpo.hdr[0])
    test = binary_KPD_model(kpo, params)
    err = kpo.kpd - test
    if kpo.kpe != None:
        err /= kpo.kpe
    return err

# =========================================================================
# =========================================================================
def binary_KPD_fit(kpo, p0):
    '''Performs a best binary fit search for the datasets.
    
    -------------------------------------------------------------
    p0 is the initial guess for the parameters 3 parameter vector
    typical example would be : [100.0, 0.0, 5.0].
    
    returns the full solution of the least square fit:
    - soluce[0] : best-fit parameters
    - soluce[1] : covariance matrix
    ------------------------------------------------------------- '''
    
    soluce = leastsq(binary_KPD_fit_residuals, 
                     p0, args=((kpo,)),
                     full_output=1)
    
    covar = soluce[1]
    return soluce

# =========================================================================
# =========================================================================
def chi2map_sep_con(kpo, pa=0, reduced=False, cmap=None,
                    srng=[10.0, 200.0, 100], crng=[10.0, 50.0, 100]):

    ''' Draws a 2D chi2 map of the binary parameter space for a given PA
    ---------------------------------------
    PA: position angle in degrees.
    other options:
    - reduced: (boolean) reduced chi2 or not
    - cmap: (matplotlib color map) "jet" or "prism" are good here
    - srng: range of separations [min, max, # steps]
    - crng: range of contrasts   [min, max, # steps]
    --------------------------------------- '''
    nullvar = np.var(kpo.kpd)
    seps = srng[0] + (srng[1]-srng[0]) * np.arange(srng[2])/float(srng[2])
    cons = crng[0] + (crng[1]-crng[0]) * np.arange(crng[2])/float(crng[2])

    chi2map = np.zeros((srng[2], crng[2]))

    for i,sep in enumerate(seps):
        for j,con in enumerate(cons):
            chi2map[i,j] = np.sum(
                binary_KPD_fit_residuals([sep, pa, con], kpo)**2)
        sys.stdout.write("\r(%3d/%d): sep = %.2f mas" % (i+1, srng[2], sep))
        sys.stdout.flush()
    if reduced:
        chi2map /= kpo.kpd.size
        nullvar /= kpo.kpd.size
    plt.clf()

    asp = 1./np.abs((srng[1]-srng[0])/(crng[1]-crng[0]))
    plt.imshow(chi2map, aspect=asp, cmap=cmap,
               extent=[crng[0], crng[1], srng[0], srng[1]])
    plt.ylabel("ang. sep. (mas)")
    plt.xlabel("contrast prim/sec")

    print("\n----------------------------------------")
    print("Best "+reduced*"reduced "+"chi2 = %.3f" % (chi2map.min(),))
    print(reduced*"reduced "+"chi2 for null hypothesis: %.3f" % (nullvar,))
    arg = chi2map.argmin()
    print("Obtained for sep = %.1f mas" % (seps[arg // chi2map.shape[1]]))
    print("Obtained for con = %.1f"     % (cons[arg  % chi2map.shape[1]]))
    print("----------------------------------------")
    return chi2map
    
# =========================================================================
# =========================================================================
def chi2map_sep_pa(kpo, con=10., reduced=False, cmap=None,
                    srng=[10.0, 200.0, 100], arng=[0.0, 360.0, 60]):

    ''' Draws a 2D chi2 map of the binary parameter space for a given contrast
    ---------------------------------------
    con: contrast (primary/secondary)
    other options:
    - reduced: (boolean) reduced chi2 or not
    - cmap: (matplotlib color map) "jet" or "prism" are good here
    - srng: range of separations     [min, max, # steps]
    - arng: range of position angles [min, max, # steps]
    --------------------------------------- '''

    nullvar = np.var(kpo.kpd)
    seps = srng[0] + (srng[1]-srng[0]) * np.arange(srng[2])/float(srng[2])
    angs = arng[0] + (arng[1]-arng[0]) * np.arange(arng[2])/float(arng[2])
    nkphi = kpo.kpi.nkphi

    chi2map = np.zeros((srng[2], arng[2]))

    for i,sep in enumerate(seps):
        for j,th in enumerate(angs):
            chi2map[i,j] = np.sum(
                binary_KPD_fit_residuals([sep, th, con], kpo)**2)
        sys.stdout.write("\r(%3d/%d): sep = %.2f mas" % (i+1, srng[2], sep))
        sys.stdout.flush()

    if reduced:
        chi2map /= kpo.kpd.size
        nullvar /= kpo.kpd.size
    plt.clf()

    asp = 1./np.abs((srng[1]-srng[0])/(arng[1]-arng[0]))

    plt.imshow((chi2map), aspect=asp, cmap=cmap,
               extent=[arng[0], arng[1], srng[0], srng[1]])
    plt.xlabel("Position Angle (deg)")
    plt.ylabel("Separation (mas)")

    print("\n----------------------------------------")
    print("Best "+reduced*"reduced "+"chi2 = %.3f" % (chi2map.min()))
    print(reduced*"reduced "+"chi2 for null hypothesis: %.3f" % (nullvar,))
    arg = chi2map.argmin()
    print("Obtained for sep = %.1f mas" % (seps[arg // chi2map.shape[1]]))
    print("Obtained for ang = %.1f deg" % (angs[arg  % chi2map.shape[1]]))
    print("----------------------------------------")
    return chi2map
    
# =========================================================================
# =========================================================================
def super_cal_coeffs(src, cal, model=None, regul="None"):
    ''' Determine the best combination of vectors in the "cal" array of
    calibrators for the source "src". 

    Regularisation is an option:
    - "None"         ->  no regularisation
    - anything else  ->  Tikhonov regularisation  '''

    A      = np.matrix(cal.kpd).T # vector base matrix
    ns     = A.shape[1]           # size of the vector base
    b      = src.kpd              # column vector
    if model != None: b -= model  # optional model subtraction

    if regul == "None":
        coeffs = np.dot(np.linalg.pinv(np.dot(A.T,A)),
                        np.dot(A.T,b).T)
    else:
        coeffs = np.dot(np.linalg.pinv(np.dot(A.T,A)+np.identity(ns)),
                        np.dot(A.T,b).T)
    return coeffs.T

# =========================================================================
# =========================================================================
def cal_search(src, cal, regul="None"):
    ''' Proceed to an exhaustive search of the parameter space.

    In each point of the space, the best combination of calibrator frames
    is found, and saved. Not entirely sure what this is going to be useful
    for...
    '''
    
    ns, s0, s1 = 50, 20.0, 200.0
    nc, c0, c1 = 50,  2.0, 100.0
    na, a0, a1 = 60,  0.0, 360.0
    
    seps = s0 + np.arange(ns) * (s1-s0) / ns
    angs = a0 + np.arange(na) * (a1-a0) / na
    cons = c0 + np.arange(nc) * (c1-c0) / nc

    coeffs = super_cal_coeffs(src, cal, None, "tik")
    cals = np.zeros((ns, na, nc, coeffs.size))

    for i,sep in enumerate(seps):
        for j, ang in enumerate(angs):
            for k, con in enumerate(cons):
                model  = binary_model([sep, ang, con], src.kpi, src.hdr)
                coeffs = super_cal_coeffs(src, cal, model, "tik")
                cals[i,j,k] = coeffs

    return cals

# =========================================================================
# =========================================================================
def chi2_volume(src, cal=None, regul="None"):
    ''' Proceed to an exhaustive search of the parameter space
    '''

    ns, s0, s1 = 50, 20.0, 200.0
    nc, c0, c1 = 50,  2.0, 100.0
    na, a0, a1 = 60,  0.0, 360.0
    
    seps = s0 + np.arange(ns) * (s1-s0) / ns
    angs = a0 + np.arange(na) * (a1-a0) / na
    cons = c0 + np.arange(nc) * (c1-c0) / nc

    chi2vol = np.zeros((ns,na,nc))

    sig = src.kpd
    if cal != None: sig -= cal.kpd

    for i,sep in enumerate(seps):
        for j, ang in enumerate(angs):
            for k, con in enumerate(cons):
                res = binary_multi_KPD_fit_residuals([sep, ang, con], src)
                chi2vol[i,j,k] = ((res)**2).sum()

    return chi2vol

# =========================================================================
# =========================================================================
def correlation_plot(kpo, params=[250., 0., 5.],b=None, plot_error=True):
    '''Correlation plot between KP object and a KP binary model
    
    Parameters are:
    --------------
    - kpo: one instance of kernel-phase object
    - params: a 3-component array describing the binary (sep, PA and contrast)

    Option:
    - plot_error: boolean, errorbar or regular plot
    --------------------------------------------------------------------------
    '''
    if b == None:

        mm = np.round(np.max(np.abs(kpo.kpd)), -1)
        
        f1 = plt.figure()
        sp0 = f1.add_subplot(111)
        if plot_error:
            sp0.errorbar(binary_KPD_model(kpo, params), kpo.kpd, 
                         yerr=kpo.kpe, linestyle='None')
        else:
            sp0.plot(binary_KPD_model(kpo, params), kpo.kpd, 'bo')
        sp0.plot([-mm,mm],[-mm,mm], 'g')
        sp0.axis([-mm,mm,-mm,mm])

        rms = np.std(binary_KPD_fit_residuals(params, kpo))
        msg  = "Model:\n sep = %6.2f mas" % (params[0],)
        msg += "\n   PA = %6.2f deg" % (params[1],)
        msg += "\n  con = %6.2f" % (params[2],)        
        msg += "\n(rms = %.2f deg)" % (rms,)
                
        plt.text(0.0*mm, -0.75*mm, msg, 
                 bbox=dict(facecolor='white'), fontsize=14)
                
        msg = "Target: %s\nTelescope: %s\nWavelength = %.2f um" % (
            kpo.kpi.name, kpo.hdr['tel'], kpo.hdr['filter']*1e6)
                
        plt.text(-0.75*mm, 0.5*mm, msg,
                  bbox=dict(facecolor='white'), fontsize=14)
        
        plt.xlabel('Data kernel-phase signal (deg)')
        plt.ylabel('Kernel-phase binary model (deg)')
        plt.draw()
        return None

    else:
        plt.clf()
        mm = np.round(np.max(np.abs(kpo.kpd)), -1)
        plt.errorbar(b,kpo.kpd,yerr=kpo.kpe,fmt='b.')
        plt.plot([-mm,mm],[-mm,mm], 'g')
        plt.axis('tight',fontsize='large')
        plt.xlabel('Model Kernel Phases',fontsize='large')
        plt.ylabel('Kernel Phase Signal', fontsize='large')
        plt.title('Kernel Phase Correlation Diagram',
                  fontsize='large')
        plt.draw()
        plt.show()


# =========================================================================
# =========================================================================

def nested(kpo, prange=[30.,250., 0.,360.,1.5,5.], npoints=100, nsteps = 3000,
           njumps=20, mode='kp'):
    '''Run a nested sampling routine with parameters [sepmin, sepmax, anglemin,
    anglemax, cmin, cmax], npoints, nsteps and njumps.
    Call with argument mode = 'vis' to use vis fitting etc.
    This is sufficiently modular that with a different loglikelihood function
    you can use it straight out of the box.

    Note: this function requires properly calibrated observables.'''

    '''--------------------------------------------------------
    Generate 100 active points from a uniform prior in angle,
    and from a Jeffreys prior in separation and contrast
    (IMPLEMENT!)
    --------------------------------------------------------'''
    nparams = len(prange)/2
    pmin = np.array([prange[j*2] for j in range(nparams)])
    pmax = np.array([prange[j*2+1] for j in range(nparams)])
    print pmin,pmax
                 
    params = np.zeros((nparams,npoints))
    print params.shape
    
    for j in range(nparams): # initialise points
             params[j,:] = (pmax[j]-pmin[j])*rand(npoints)+pmin[j] 

    L = np.zeros(npoints)

    rejected = np.zeros((nsteps,nparams))
    rejL = np.zeros(nsteps)

    #lognull = -np.sum(((kpo.kpd)/a.kpe)**2)/2 +logpre

    progress = 99 # counter for displays

    for k in range(0,npoints): # calculate likelihood

        L[k] = loglikelihood(kpo,params[:,k],mode=mode)

    for j in range(0,nsteps): # loop for iterations

        print 'Commencing iteration',j

        # store rejected binary parameters and re-sort L
        jump = (pmax-pmin)/float(npoints)

        for m in range(0,nparams):
            rejected[j,m] = params[m,L==L.min()]
        
        rejL[j] = L[L==L.min()]
        
        args = np.argsort(L)
        params, L = params[:,args], L[args]

        acceptances = 0 # initialise for coming loop
        rejections = 0
        l = 0
        trys = np.zeros(nparams)

        choose = int(np.ceil(rand()*(npoints-1))) # make a random choice

        starts = params[:,choose]
        news = starts
        newl = loglikelihood(kpo,news,mode=mode)
        
        while l<njumps: # random walk to generate a new active point

            trys = news + jump*randn(nparams)
            trys[1] = np.mod(trys[1],360) # assumes second parameter is an angle
            
            like = loglikelihood(kpo,trys,mode=mode)
            
            proc = np.all(pmin<trys) and np.all(trys<pmax) and like>L[0]
            
            if proc : #accept and move on
                acceptances += 1
                news = trys
                newl = like
            else: # dwell on it
                rejections += 1

            #adjust jumps - factors from John Skilling's example

            if acceptances > rejections: # jump further if you can afford to
                jump *= np.exp(1/float(acceptances))
                
            else: # don't jump as far if you're making it worse too often
                jump /= np.exp(1/float(rejections))
                
            l += 1

            if (l >= njumps and acceptances == 0): # start somewhere else
                l=0
                choose = int(np.ceil(rand()*(npoints-1)))
                news = params[:,choose]

        params[:,0] = news
        L[0] = newl

        print l, 'steps taken', acceptances/float(l),'acceptance rate'
        print 'Sep',params[0,0],'Angle',params[1,0],'Contrast',params[2:,0]

        if progress == 99: # make a plot every 100 iterations
            print 'Generating plots'

            plt.figure(0)
            plt.clf()
            plt.scatter(params[0],params[1])
            plt.axis([pmin[0],pmax[0],pmin[1],pmax[1]])
            plt.xlabel('Separation (mas)')
            plt.ylabel('Position Angle (degrees)')
            plt.title('Active Points')
            plt.draw()
            
            #plt.show()

            plt.figure(1)
            plt.clf()
            plt.scatter(params[0],params[2])
            plt.axis([pmin[0],pmax[0],pmin[2],pmax[2]])
            plt.xlabel('Separation (mas)')
            plt.ylabel('Contrast Ratio')
            plt.title('Active Points')
            plt.draw()

            if nparams > 3:
                plt.figure(2)
                plt.clf()
                plt.scatter(params[0],params[3])
                plt.axis([pmin[0],pmax[0],pmin[3],pmax[3]])
                plt.xlabel('Separation (mas)')
                plt.ylabel('Contrast Ratio 2')
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
    logpre = 0 #workaround for the moment

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

    means = [np.sum(pdf*rejected[:,k]) for k in range(nparams)]
    ds = [np.sqrt(np.sum(pdf*(rejected[:,k] - means[k])**2)) for k in range(nparams)]

    output = means+ds+[logE,dlogE]
    print output
    
    return output
# =========================================================================
# =========================================================================

def loglikelihood(kpo,params,mode='kp'):
    '''Define the chi-squared function - omits logpre!
    Takes a structure a and parameters sep, angle and contrast.
    Pass the mode 'vis' to use visibility fitting

    If fitting over multiple kpo bands, params are listed as [sep, angle, c1,
    c2, ... cn]. '''

    if type(kpo) != list:
        [sep,angle,contrast] = params

        if mode == 'kp':
            modl_ker = binary_model([sep,angle,contrast],kpo.kpi,kpo.hdr)

            chisquared = np.sum(((modl_ker - kpo.kpd)/(kpo.kpe))**2)
        elif mode == 'vis':
            test = vis2_binary(kpo.kpi.uv[:,0],kpo.kpi.uv[:,1], kpo.hdr['filter'],
                               [sep,angle+90.0-kpo.hdr['orient'],contrast])
            chisquared = np.sum(((test - kpo.vis2)/(kpo.vis2e))**2)
        else:
            print 'Invalid mode'

    else:
        nsets = len(kpo)
        chisquareds = np.zeros(nsets)
        [sep, angle] = [params[0],params[1]]
        
        for j in range(nsets):
            contrast = params[j+2]
            
            if mode == 'kp':
                modl_ker = binary_model([sep,angle,contrast],kpo[j].kpi,kpo[j].hdr)

                chisquareds[j] = np.sum(((modl_ker - kpo[j].kpd)/(kpo[j].kpe))**2)
            elif mode == 'vis':
                test = vis2_binary(kpo[j].kpi.uv[:,0],kpo[j].kpi.uv[:,1], kpo[j].hdr['filter'],
                                   [sep,angle+90.0-kpo[j].hdr['orient'],contrast])
                chisquareds[j] = np.sum(((test - kpo[j].vis2)/(kpo[j].vis2e))**2)
            else:
                print 'Invalid mode'
        chisquared = np.sum(chisquareds)

    #logpre = a.nkphi*np.sum(-1/2 * np.log(2*np.pi) - np.log(a.kpe))
    #logpre = 0 
    #like = logpre-chisquared/2
    like = -chisquared/2
    
    return like
