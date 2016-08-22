''' --------------------------------------------------------------------
    pyker is a python module to create, and extract Kernel-phase data 
    structures, using the theory of Martinache, 2010, ApJ, 724, 464.
    ----

    The module is constructed around two main classes:
    -------------------------------------------------

    - KernelPhase: object that contains the linear model for the 
      optical system of interest. Properties of this model are:
      --> name   : name of the model (HST, Keck, Annulus_19, ...)
      --> mask   : array of coordinates for pupil sample points
      --> uv     : matching array of coordinates in uv plane (baselines)
      --> RED    : vector coding the redundancy of these baselines
      --> TFM    : transfer matrix, linking pupil-phase to uv-phase
      --> KerPhi : array storing the kernel-phase relations

    - KerPhase_Data: objects that contains Ker-phase data extracted
      from actual images, using the model of KerPhase_Relation + some
      other information: plate scale, wavelength, ...
      -------------------------------------------------------------------- '''

#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt, pyfits as pf
import pickle, os, pdb, matplotlib.cm as cm
import pylab
from scipy.special import j1
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

from pyker_tools import *

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

dtor = np.pi/180.0

# -------------------------------------------------
# set some defaults to display images that will
# look more like the DS9 display....
# -------------------------------------------------
#plt.set_cmap(cm.gray)
#(plt.rcParams)['image.origin']        = 'lower'
#(plt.rcParams)['image.interpolation'] = 'nearest'
# -------------------------------------------------

class KernelPhase(object):
    ''' Fundamental kernel-phase relations

    -----------------------------------------------------------------------
    This object condenses all the knowledge about a given instrument pupil 
    geometry into a series of arrays useful for kernel-phase analysis as 
    well as for other purposes, such as wavefront sensing.
    ----------------------------------------------------------------------- '''

    name = "" # default array name. Should be descriptive of the array geometry

    # =========================================================================
    # =========================================================================

    def __init__(self, file=None):
        ''' Default instantiation of a KerPhase_Relation object:

        -------------------------------------------------------------------
        Default instantiation of this KerPhase_Relation class is achieved
        by loading a pre-made file, containing all the relevant information
        -------------------------------------------------------------------'''
        try:
            # -------------------------------
            # load the pickled data structure
            # -------------------------------
            myf = open(file, "r")
            data = pickle.load(myf)
            myf.close()

            # -------------------------------
            # restore the variables for this 
            # session of Ker-phase use!
            # -------------------------------
            #print "TEST:", data['toto']

            try:    self.name = data['name']
            except: self.name = "UNKNOWN"

            try:    self.file = file
            except: self.file = 'File not specified'

            self.uv     = data['uv']
            self.mask   = data['mask']
            self.RED    = data['RED']
            self.KerPhi = data['KerPhi']
            self.TFM    = data['TFM']
        
            self.nbh   = self.mask.shape[0]
            self.nbuv  = self.uv.shape[0]
            self.nkphi = self.KerPhi.shape[0]
        
        except: 
            print("File %s isn't a valid Ker-phase data structure" % (file))
            try: self.from_coord_file(file)
            except:
                print("Failed.")
                return None

    # =========================================================================
    # =========================================================================

    def from_coord_file(self, file, array_name=""):
        ''' Creation of the KerPhase_Relation object from a pupil mask file:

        ----------------------------------------------------------------
        This is the core function of this class, really...

        Input is a pupil coordinates file, containing one set of (x,y) 
        coordinates per line. Coordinates are in meters. From this, all 
        the intermediate products that lead to the kernel-phase matrix 
        KerPhi are calculated.
        ---------------------------------------------------------------- '''
        self.mask = 1.0 * np.loadtxt(file) # sub-Ap. coordinate files 
        self.nbh  = self.mask.shape[0]   # number of sub-Ap

        ndgt = 6 # number of digits of precision for rounding
        prec = 10**(-ndgt)

        # ================================================
        # Determine all the baselines in the array.
        # ================================================

        # 1. Start by doing all the possible combinations of coordinates 
        # --------------------------------------------------------------
        # in the array to calculate the baselines. The intent here, is 
        # to work with redundant arrays of course, so there will be plenty 
        # of duplicates.

        nbh = self.nbh # local representation of the class variable
        uvx = np.zeros(nbh * (nbh-1)) # prepare empty arrays to store
        uvy = np.zeros(nbh * (nbh-1)) # the baselines

        k = 0 # index for possible combinations (k = f(i,j))
        
        uvi = np.zeros(nbh * (nbh-1), dtype=int) # arrays to store the possible
        uvj = np.zeros(nbh * (nbh-1), dtype=int) # combinations k=f(i,j) !!


        for i in range(nbh):     # do all the possible combinations of
            for j in range(nbh): # sub-apertures
                if i != j:
                    uvx[k] = self.mask[i,0] - self.mask[j,0]
                    uvy[k] = self.mask[i,1] - self.mask[j,1]
                    # ---
                    uvi[k], uvj[k] = i, j
                    k+=1

        a = np.unique(np.round(uvx, ndgt)) # distinct u-component of baselines
        nbx    = a.shape[0]                # number of distinct u-components
        uv_sel = np.zeros((0,2))           # array for "selected" baselines

        for i in range(nbx):     # identify distinct v-coords and fill uv_sel
            b = np.where(np.abs(uvx - a[i]) <= prec)
            c = np.unique(np.round(uvy[b], ndgt))
            nby = np.shape(c)[0] # number of distinct v-compoments
            for j in range(nby):
                uv_sel = np.append(uv_sel, [[a[i],c[j]]], axis=0)

        self.nbuv = np.shape(uv_sel)[0]/2 # actual number of distinct uv points
        self.uv   = uv_sel[:self.nbuv,:]  # discard second half (symmetric)
        print "%d distinct baselines were identified" % (self.nbuv,)

        # 2. Calculate the transfer matrix and the redundancy vector
        # --------------------------------------------------------------
        self.TFM = np.zeros((self.nbuv, self.nbh), dtype=float) # matrix
        self.RED = np.zeros(self.nbuv, dtype=float)             # Redundancy

        for i in range(self.nbuv):
            a=np.where((np.abs(self.uv[i,0]-uvx) <= prec) *
                       (np.abs(self.uv[i,1]-uvy) <= prec))
            self.TFM[i, uvi[a]] +=  1.0
            self.TFM[i, uvj[a]] += -1.0
            self.RED[i]         = np.size(a)

        # 3. Determine the kernel-phase relations
        # ----------------------------------------

        # One sub-aperture is taken as reference: the corresponding
        # column of the transfer matrix is discarded. TFM is now a
        # (nbuv) x (nbh - 1) array.
        
        # The choice is up to the user... but the simplest is to
        # discard the first column, that is, use the first aperture
        # as a reference?

        self.TFM = self.TFM[:,1:] # cf. explanation
        U, S, Vh = np.linalg.svd(self.TFM.T, full_matrices=1)

        S1 = np.zeros(self.nbuv)
        S1[0:nbh-1] = S

        self.nKPhi  = np.size(np.where(abs(S1) < 1e-3)) # number of Ker-phases
        KPhiCol     = np.where(abs(S1) < 1e-3)[0]
        self.KerPhi = np.zeros((self.nKPhi, self.nbuv)) # allocate the array

        for i in range(self.nKPhi):
            self.KerPhi[i,:] = (Vh)[KPhiCol[i],:]

        print '-------------------------------'
        print 'Singular values for this array:\n', np.round(S, ndgt)
        print '\nRedundancy Vector:\n', self.RED
        self.name = array_name

    # =========================================================================
    # =========================================================================

    def plot_pupil_and_uv(self, xymax = 8.0):
        ''' Nice plot of the pupil sampling and matching uv plane.

        --------------------------------------------------------------------
        xymax just specifies the size of the region represented in the plot,
        expressed in meters. Should typically be slightly larger than the 
        largest baseline in the array.
        --------------------------------------------------------------------'''

        plt.clf()
        f0 = plt.subplot(121)
        f0.plot(self.mask[:,0], self.mask[:,1], 'bo')
        f0.axis([-xymax, xymax, -xymax, xymax])
        plt.title(self.name+' pupil')
        f1 = plt.subplot(122)

        f1.plot(self.uv[:,0],   self.uv[:,1], 'bo') # plot baselines + symetric
        f1.plot(-self.uv[:,0], -self.uv[:,1], 'ro') # for a "complete" feel
        plt.title(self.name+' uv coverage')
        f1.axis([-2*xymax, 2*xymax, -2*xymax, 2*xymax])

        # complete previous plot with redundancy of the baseline
        # -------------------------------------------------------
        dy = 0.1*abs(self.uv[0,1]-self.uv[1,1]) # to offset text in the plot.
        for i in range(self.nbuv):
            f1.text(self.uv[i,0]+dy, self.uv[i,1]+dy, 
                    int(self.RED[i]), ha='center')

        plt.draw()

    # =========================================================================
    # =========================================================================

    def save_to_file(self, file):
        ''' Export the KerPhase_Relation data structure into a pickle
        
        ----------------------------------------------------------------
        Nothing really noteworthy here: the data structure is simply
        saved for future use.
        ----------------------------------------------------------------  '''
        try: data = {'name'   : self.name,
                     'mask'   : self.mask,
                     'uv'     : self.uv,
                     'TFM'    : self.TFM,
                     'KerPhi' : self.KerPhi,
                     'RED'    : self.RED}
        except:
            print("KerPhase_Relation data structure is incomplete")
            print("File %s wasn't saved!" % (file,))
            return None
        # -------------
        try: myf = open(file, "w")
        except:
            print("File %s cannot be created."+
                  " KerPhase_Relation data structure wasn't saved." % (file,))
            return None
        # -------------
        pickle.dump(data, myf)
        myf.close()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KernelPhaseData(KernelPhase):
    ''' Class used to manipulate Kernel-phase data (actual or modeled)

        ------------------------------------------------------------------
        Class inherits from KernelPhase.
        ------------------------------------------------------------------ '''

    def __init__(self, KerPhiFname):
        '''Default instantiation from a Kernel-phase data file.

        ------------------------------------------------------------
        Uses parent class KernelPhase init for KerPhase information,
        completes it with KerPhaseData specifics
        ------------------------------------------------------------ '''
        super(KernelPhaseData, self).__init__(KerPhiFname)
        self.load_KerPhase_data(KerPhiFname)

    # =========================================================================
    # =========================================================================

    def load_KerPhase_data(self, KerPhiFname):
        '''Load the KerPhaseData specific information.

        If the file happens to be incomplete, the object will contain 
        a bunch of empty arrays.'''
        
        try:
            myf = open(KerPhiFname, "r")
            data = pickle.load(myf)
            myf.close()
        except:
            print("File %s isn't a valid Ker-phase data structure" % (file))
            return None

        try:    self.comments = data['comments']
        except: self.comments = ""

        try:    self.filter = data['filter']
        except: self.filter = 0.0

        try:    self.kp_signal = data['kp_signal']
        except: self.kp_signal = np.zeros(self.KerPhi.shape[0])

        try:    self.kp_error = data['kp_error']
        except: self.kp_error = np.zeros(self.KerPhi.shape[0])

        
    # =========================================================================
    # =========================================================================

    def binary_fit(self, p0):
        '''Performs a best binary fit search for the dataset.

        -------------------------------------------------------------
        p0 is the initial guess for the parameters 3 parameter vector
        typical example would be : [100.0, 0.0, 5.0].

        returns the full solution of the least square fit:
        - soluce[0] : best-fit parameters
        - soluce[1] : covariance matrix
        ------------------------------------------------------------- '''

        if np.all(self.kp_signal == 0.0):
            print("KernelPhaseData instance contains no data.\n" +
                  "Model fit is not possible\n")
            return None
        
        if np.all(self.kp_error == 0.0):
            print("KernelPhaseData instance is not calibrated.\n")
            soluce = leastsq(self.bin_fit_residuals, p0, args=(self.kp_signal),
                         full_output=1)
        else: soluce = curve_fit(self.kphitest,range(0,self.nkphi),self.kp_signal,p0,sigma=self.kp_error)

        self.covar = soluce[1]
        soluce[0][1] = np.mod(soluce[0][1],360.) # to get consistent position angle measurements
        return soluce

    # =========================================================================
    # =========================================================================

    def bin_fit_residuals(self, params, data):
        p = np.copy(params)
        p[1] += 90.0 - self.info['orient']
        filter = self.info['filter']

        testPhi = phase_binary(self.uv[:,0], self.uv[:,1], filter, p)
        test = np.dot(self.KerPhi, self.RED * testPhi)

        err = (data - test)
        return err

    # =========================================================================
    # =========================================================================

    def binary_model(self, params, kp_info):
        ''' Creates a binary Kernel-phase model.
        
        ------------------------------------------------------------------ 
        uses a simple 5 parameter binary star model for the uv phases that
        should be observed with the provided geometry.

        The additional parameter is a reference to a KernelPhaseData info
        structure, that should for now, contain valid information for:
        - filter (wavelength of the filter in meters)
        - orient (frame orientation)
        ------------------------------------------------------------------ '''

        params2 = np.copy(params)
        params2[1] += 90.0 - kp_info['orient'] # 
        filter = kp_info['filter']

        testPhi = phase_binary(self.uv[:,0], self.uv[:,1], filter, params2)


        self.kp_signal = np.dot(self.KerPhi, self.RED * testPhi)
        self.kp_error  = np.zeros(self.KerPhi.shape[0])
        self.comments += "Binary model, parameters = %s\n" % (str(params))

        # -----------------------------------------------------------------
        # fill up the info structure with bogus numbers (except for filter)
        # -----------------------------------------------------------------
        data = {
            'tel'    : 'PERFECT',
            'pscale' : 0.0,
            'fname'  : 'SIMULATION',
            'odate'  : 'today',  # UTC date of observation
            'otime'  : 'now',    # UTC time of observation
            'tint'   : 0.0,      # integration time (sec)
            'coadds' : 1,        # number of coadds
            'RA'     : 0.0,      # right ascension (deg)
            'DEC'    : 0.0,      # declination (deg)
            'filter' : filter,   # central wavelength (meters)
            'orient' : 0.0,      # position angle of the frame (deg)
            'params' : params,   # parameters of the simulated binary
            'covar'  : np.zeros((3,3)) # covariance matrix for fit
            }
        self.info = data
        return self.kp_signal

    # =========================================================================
    # =========================================================================

    def kphitest(self,index,sep,angle,contrast):
        ''' Dummy taken out of binary_model to play with in error-weighted fitting'''
        params = [sep,angle,contrast]
        params2 = np.copy(params)
        params2[1] += 90.0 - self.info['orient'] # 
        filter = self.info['filter']

        testPhi = phase_binary(self.uv[:,0], self.uv[:,1], filter, params2)

        kphasemodel = np.dot(self.KerPhi, self.RED * testPhi)

        return kphasemodel[index]

    # =========================================================================
    # =========================================================================

    def fit_map(self,sep=None,angle=None,contrast=None,
                sepmin = 40.0,sepmax = 1000.0,anglemin=0.0, anglemax=359.9,
                cmin = 1.01,cmax = 10.0, random = None, display= True,
                nclevs = 20,clevsmax=12,filename=None,returnmap=False):
        ''' Takes two or more of the binary parameters, calculates a range in
        the other parameters, and plots a chi-squared map over this space. Include
        random != None for a random grid rather than a regular grid.

        This starts out with a projection. To go over all three parameters,
        you should use chi2_3d'''

        filter = self.info['filter']

        # first make the parameter space

        nx,ny = 120,120

        if random is not None: # calculate a random grid
            if sep is not None:
                key = 1
                ttl = 'Reduced chi-squared map, sep = ' + str(sep)
                xaxis = anglemin + np.sort(np.random.rand(nx))*(anglemax-anglemin)
                xlbl = 'Position Angle (degrees)'
                xmin, xmax = anglemin, anglemax
                yaxis = cmin + np.sort(np.random.rand(ny))*(cmax-cmin)
                ylbl = 'Contrast'
                ymin, ymax = cmin, cmax
                xy = np.meshgrid(xaxis,yaxis)
            elif angle is not None:
                key = 2
                ttl = 'Reduced chi2 map, position angle = '+str(angle)
                xaxis = sepmin + np.sort(np.random.rand(nx))*(sepmax-sepmin)
                xlbl = 'Separation (mas)'
                xmin, xmax = sepmin,sepmax
                yaxis = cmin + np.sort(np.random.rand(ny))*(cmax-cmin)
                ylbl = 'Contrast'
                ymin, ymax = cmin, cmax
            elif contrast is not None:
                key = 3
                ttl = 'Reduced chi-squared map, contrast = ' + str(contrast)
                xaxis = anglemin + np.sort(np.random.rand(nx))*(anglemax-anglemin)
                xlbl = 'Position Angle (degrees)'
                xmin, xmax = anglemin, anglemax
                yaxis = sepmin + np.sort(np.random.rand(ny))*(sepmax-sepmin)
                ylbl = 'Separation (mas)'
                ymin, ymax = sepmin, sepmax
            else:
                print "fit_map passed invalid argument(s)"
        else: #calculate a regular grid
            if sep is not None:
                key = 1
                ttl = 'Reduced chi-squared map, sep = ' + str(sep)
                xaxis = np.linspace(anglemin,anglemax,num=nx)
                xlbl = 'Position Angle (degrees)'
                xmin, xmax = anglemin, anglemax
                yaxis = np.linspace(cmin,cmax,num=ny)
                ylbl = 'Contrast'
                ymin, ymax = cmin, cmax
            elif angle is not None:
                key = 2
                ttl = 'Reduced chi2 map, position angle = '+str(angle)
                angle = angle
                xaxis = np.linspace(sepmin,sepmax,num=nx)
                xlbl = 'Separation (mas)'
                xmin, xmax = sepmin,sepmax
                yaxis = np.linspace(cmin,cmax,num=ny)
                ylbl = 'Contrast'
                ymin, ymax = cmin, cmax
            elif contrast is not None:
                key = 3
                ttl = 'Reduced chi2 map, contrast = ' + str(contrast)
                xaxis = np.linspace(anglemin,anglemax,num = nx)
                xlbl = 'Position Angle (degrees)'
                xmin, xmax = anglemin, anglemax
                yaxis = np.linspace(sepmin,sepmax,num=ny)
                ylbl = 'Separation (mas)'
                ymin, ymax = sepmin, sepmax
            else:
                print "fit_map passed invalid argument(s)"

        # generate a uv phase map for each of the parameter sets

        xx,yy = np.meshgrid(xaxis,yaxis)

        #calculate models

        if key ==1:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[sep,xx[j][i]+90.0 -self.info['orient'],yy[j][i]])
                       for i in range(0,nx)] for j in range(0,ny)]
        elif key ==2:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[xx[j][i],angle+90.0 -self.info['orient'],yy[j][i]])
                       for i in range(0,nx)] for j in range(0,ny)]
        elif key ==3:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[yy[j][i],xx[j][i]+90.0 -self.info['orient'],contrast])
                       for i in range(0,nx)] for j in range(0,ny)]
            
        #calculate kerphases
            
        kphasemodels = [[np.dot(self.KerPhi,self.RED*models[j][i]) for i in range (0,nx)]
                        for j in range(0,ny)]

        #calculate departures from the model and the null hypothesis

        diff = [[np.divide(np.subtract(self.kp_signal,kphasemodels[j][i])**2,
                           self.kp_error**2) for i in range(0,nx)] for j in range(0,ny)]

        null = np.divide(self.kp_signal**2,self.kp_error**2)
        
        #calculate chi2 for the model and the null hypothesis

        chi2 = [[np.sum(diff[j][i]) for i in range(0,nx)] for j in range(0,ny)]
        
        chi2null= np.sum(null)

        # calculate reduced chi2

        dof = float(self.nkphi - 4)

        chi2 = np.array(chi2)/dof

        chi2null = np.array(chi2null)/dof

        # calculate significance

        bestfit = chi2.min()

        significance = np.sqrt(chi2null*dof - bestfit*dof)

        # find the best-fit coordinates

        bestindices = np.unravel_index(chi2.argmin(),chi2.shape)

        bestx, besty = xx[bestindices], yy[bestindices]

        if display is True: #default, plot it
            print 'Best chi-squared = ', bestfit

            print 'Null hypothesis chi-squared = ', chi2null

            print 'Significance estimated at ', significance

            print 'Best fit coordinates = ', bestx, besty

            #create plot

            clevs = chi2.min()*np.linspace(1,clevsmax,nclevs)
            plt.clf()

            plt.contourf(xx,yy,chi2,levels = clevs,cmap=plt.cm.bone,
                         antialiased = True)
            plt.axis('tight')
            plt.xlabel(xlbl)
            plt.ylabel(ylbl)
            plt.title(ttl)
            plt.colorbar()
            
            if filename is not None:
                plt.savefig(filename+'.png',format='png')
            else:
                plt.show()
            
            return None
        
        elif returnmap is False: # return parameters in a string p = [bestfit,significance,bestx,besty]
            p = [bestfit,significance,bestx,besty]
            return p

        elif returnmap is True: # return the chi2 map itself
            return chi2*dof, xx, yy

    # =========================================================================
    # =========================================================================

    def adaptive_fit(self,contrast,sepmin,sepmax,tolerance=1.0):
        '''Do an adaptive grid fit with a given chi-squared tolerance, a guess
        for the contrast seed, and a guess for the separation range.'''
        count = 0
        chi2 = tolerance*5 # dummy initial value
        sepdiff = sepmax-sepmin
        anglemin = 0
        anglemax = 360
        anglediff = anglemax-anglemin

        print "About to commence adaptive grid fit"

        for j in range(0,6):
            p  = self.fit_map(contrast=contrast,sepmin=sepmin,sepmax=sepmax,
                           anglemin=anglemin,anglemax=anglemax,display=False)
            #sepmin = p[3]-sepdiff/2.
            #sepmax = p[3]+sepdiff/2.
            #sepdiff = sepmax-sepmin
            anglemin = p[2]-anglediff/2.
            anglemax = p[2]+anglediff/2.
            anglediff= anglediff/2
            chi2 = p[0]
            count+=1
            print "Separation = %.2f , Position angle = %.1f, Chi-squared = %.3f" % (p[3],p[2],chi2)

        print "Position-angle fit completed"
        print "Chi-squared =",chi2
        
        angle = p[2]          

        cmin = 1.1
        cmax = contrast*2
        cdiff = cmax-cmin
        
        for j in range(0,6):
            p = self.fit_map(angle=angle,sepmin=sepmin,sepmax=sepmax,
                             cmin=cmin,cmax=cmax,display=False)
            sepmin = p[2]-sepdiff/2.
            sepmax = p[2]+sepdiff/2.
            sepdiff = sepdiff/2
            cmin = p[3]-cdiff/2.
            cmax = p[3] + cdiff/2.
            cdiff = cdiff/2.
            chi2 = p[0]
            count+=1
            print "Separation = %.3f, Contrast = %.3f, Chi-squared = %.3f" % (p[2],p[3],chi2)
                
        print "Separation-Contrast fit completed"
        print "Chi-squared =", chi2
    
    # =========================================================================
    # =========================================================================

    def random_fit(self,sepmin,sepmax,anglemin,anglemax,cmin,cmax,n,binfit=None):
        '''Generates n random linear regression fits within the space defined by
        the arguments, then returns the chi-squared of the model, the null chi-
        squared, significance, and best-fit parameters.

        NOTE: Because the binary_fit method has a habit of running away around
        parameter space, this is at present a pretty bad function for fitting. '''

        from numpy.random import rand

        filter = self.info['filter']
        dof = float(self.nkphi - 4)

        #initialise arrays to hold the solutions and chi-squared values

        sols = list()
        chi2 = np.empty(n)

        #find the null hypothesis chi-squared

        null = np.divide(self.kp_signal**2,self.kp_error**2)
        chi2null = np.sum(null)/dof

        if binfit is not None: # if you flagged it to use the binary_fit method
            print 'Initiating binary_fit sequence'
            for j in range(0,n):
                sep = sepmin+rand()*(sepmax-sepmin)
                angle = anglemin+rand()*(anglemax-anglemin)
                contrast = cmin+rand()*(anglemax-anglemin)
                params = [sep,angle+90.0-self.info["orient"],contrast]
                try:
                    fit,covar = self.binary_fit(params)
                    testPhi = phase_binary(self.uv[:,0], self.uv[:,1], filter, fit)
                    test_signal = np.dot(self.KerPhi, self.RED * testPhi)
                    chi2[j] = np.sum(np.divide(np.subtract(self.kp_signal,test_signal)**2,
                                               self.kp_error**2))/dof
                    sols.append(fit)
                except:
                    chi2[j] = 1000
                    sols.append([0,0,0])
                    print 'Fit', j,'failed to converge'
        else: # just take the raw parameters and don't try to optimise
            for j in range(0,n):
                sep = sepmin+rand()*(sepmax-sepmin)
                angle = anglemin+rand()*(anglemax-anglemin)
                contrast = cmin+rand()*(anglemax-anglemin)
                params = [sep,angle+90.0-self.info["orient"],contrast]
                testPhi = phase_binary(self.uv[:,0], self.uv[:,1], filter, params)
                test_signal = np.dot(self.KerPhi, self.RED * testPhi)
                chi2[j] = np.sum(np.divide(np.subtract(self.kp_signal,test_signal)**2,
                                           self.kp_error**2))/dof
                sols.append(params)


        bestchi2 = chi2.min()
        bestparams = sols[chi2.argmin()]
        bestparams[1] = bestparams[1]-90.0+self.info["orient"]
        significance = np.sqrt(chi2null-bestchi2)

        print "Best chi-squared = ", bestchi2
        print "Null hypothesis chi-squared = ", chi2null
        print "Significance = ", significance
        print "Best fit parameters = ", bestparams

        return None
        
    # =========================================================================
    # =========================================================================

    def chi2_3d(self,smin,smax,anglemin,anglemax,cmin,cmax,
                   nsep=100,nangle=360,nc=100):

        ''' Generates a chi2 grid in all three variables. Cannibalised from
        Frantz' code to do this in pyapm/plot_search. Goes over three variables,
        so takes longer than fit_map, which is just for looking at projections
        quickly.

        The indices are [i,j,k] over [sep,angle,contrast].'''

        filter = self.info['filter']
        
        dof = float(self.kp_signal.shape[0] - 3)

        chi2 = np.zeros((nsep, nangle, nc))
        
        nsep = float(nsep) # make sure you're gridding with floats!
        nangle = float(nangle)
        nc = float(nc)

        # define our grid points

        seps = smin + (smax-smin) * np.arange(nsep)/nsep
        angles  = anglemin + (anglemax-anglemin) * np.arange(nangle)/nangle
        cons = cmin + (cmax-cmin) * np.arange(nc)/nc

        # Calculate null hypothesis chi-squared

        nullchi2 = np.sum((self.kp_signal/self.kp_error)**2)

        for i,sep in enumerate(seps):
            print i, sep
            for j,angle in enumerate(angles):
                for k,c in enumerate(cons):
                    test = phase_binary(self.uv[:,0], self.uv[:,1],
                                        filter, [sep,angle+90.0-self.info['orient'],c])
                    modl_ker = np.dot(self.KerPhi, self.RED*test)
                    chi2[i,j,k] = np.sum(((modl_ker-self.kp_signal)/self.kp_error)**2)

        #wherevar = np.where(vari == vari.min())

        wh = np.where(chi2 == chi2.min())

        bestsep, bestangle, bestcon = seps[wh[0][0]], angles[wh[1][0]], cons[wh[2][0]]

        #chi2 = chi2*float(dof)/chi2.min()

        #self.kp_error = self.kp_error*np.sqrt(chi2.min()/float(dof))

        #print 'Best-fit separation = %.2f, position angle = %.2f, contrast = %.2f' % (bestsep,bestangle,bestcon)

        #print 'Best chi2 obtained = %.3f' % chi2.min()

        return chi2,nullchi2

    # =========================================================================
    # =========================================================================

    def marginalplot(self,axis,smin,smax,anglemin,anglemax,cmin,cmax,
                     nsep=100,nangle=360,nc=100,chi2=None,save=None):

        ''' Plot a likelihood function marginalised over an axis.

        0 = separation
        1 = position angle
        2 = contrast ratio '''

        if chi2 is None:
            chi2 = self.chi2_3d(smin,smax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc)

        likelihood = np.exp(-chi2/2.0)
        likelihood = likelihood/(likelihood.sum())

        marg = np.sum(likelihood,axis=axis)
        marg = marg/marg.sum()

        if axis == 0:
            xx = anglemin + (anglemax-anglemin) * np.arange(nangle)/nangle
            yy = cmin + (cmax-cmin) * np.arange(nc)/nc
            xlbl = 'Position angle (deg)'
            ylbl = 'Contrast Ratio'
            ttl = 'Likelihood function, marginalised over separation'
            marg = marg.transpose() # to get the shape right!
        elif axis == 1:
            xx = smin + (smax-smin) * np.arange(nsep)/nsep
            yy = cmin + (cmax-cmin) * np.arange(nc)/nc
            xlbl = 'Separation (mas)'
            ylbl = 'Contrast Ratio'
            ttl = 'Likelihood function, marginalised over angle'
            marg = marg.transpose() # to get the shape right!
        elif axis == 2:
            xx = anglemin + (anglemax-anglemin) * np.arange(nangle)/nangle
            yy = smin + (smax-smin) * np.arange(nsep)/nsep
            xlbl = 'Position Angle (deg)'
            ylbl = 'Separation (mas)'
            ttl = 'Likelihood function, marginalised over contrast'
        else: print 'axis must be 0, 1 or 2'

        plt.pcolor(xx,yy,-marg,cmap=plt.cm.bone,
                     antialiased = True)
        
        plt.axis('tight')
        plt.xlabel(xlbl)
        plt.ylabel(ylbl)
        plt.title(ttl)
        plt.draw()
        if save is None:
            plt.show()
        elif save is not None:
            plt.savefig(save,format='png')
        else: print "Failed to specify a good filename!"

        
    # =========================================================================
    # =========================================================================


    def correlation_plot(self, model):
        '''Correlation plot between Ker-phase data and provided model'''
        mm = np.round(np.max(np.abs(self.kp_signal)), -1)

        f1 = plt.figure()
        sp0 = f1.add_subplot(111)
        sp0.errorbar(model.kp_signal, self.kp_signal,yerr=self.kp_error, fmt='b.')
        sp0.plot([-mm,mm],[-mm,mm], 'g')
        #sp0.axis([-150,150,-150,150])
        sp0.axis('tight')
        try:do_it = True - np.all(self.covar == 0.0)
        except:
            pass
        try: 
            p = model.info['params']
            self.rms = np.std(self.bin_fit_residuals(p, self.kp_signal))
            if do_it:
                msg  = "Model:\n sep=%6.2f +/- %5.2f mas" % \
                    (p[0], self.covar[0,0]*self.rms**2)
                msg += "\n  PA=%6.2f +/- %5.2f deg" % \
                    (p[1], self.covar[1,1]*self.rms**2)
                msg += "\n    c=%6.2f +/- %5.2f" % \
                    (p[2], self.covar[2,2]*self.rms**2)
            else:
                msg  = "Model:\n sep=%6.2f" % (p[0])
                msg += "\n   PA=%6.2f" % (p[1])
                msg += "\n  c=%6.2f" % (p[2])

            msg += "\n(rms= %.2f deg)" % (self.rms,)

            plt.text(0.0*mm, -0.75*mm, msg, 
                      bbox=dict(facecolor='white'), fontsize=14)

            msg = "Target: %s\nTelescope: %s\nWavelength = %.2f um" % (
                self.info['fname'], self.info['tel'], self.info['filter']*1e6)

            plt.text(-0.75*mm, 0.5*mm, msg,
                      bbox=dict(facecolor='white'), fontsize=14)
        except:
            pass

        plt.xlabel('Model kernel-phase signal (deg)')
        plt.ylabel('Kernel-phase data (deg)')
        plt.draw()
        plt.show()

    # =========================================================================
    # =========================================================================

    def residual_plot(self,model):
        '''Plot the residuals between ker-phase data and a model. Code cannibal-
        ised from correlation_plot. '''

        mm = np.round(np.max(np.abs(self.kp_signal)), -1)
        nn = np.round(np.max(np.abs(np.subtract(self.kp_signal,model.kp_signal))), -1)

        f1 = plt.figure()
        sp0 = f1.add_subplot(111)
        sp0.plot(self.kp_signal, np.subtract(self.kp_signal,model.kp_signal), 'b.')
        sp0.plot([-mm,mm],[0,0], 'g')
        sp0.axis([-mm,mm,-nn,nn])

        do_it = True - np.all(self.covar == 0.0)
        try: 
            p = model.info['params']
            self.rms = np.std(self.bin_fit_residuals(p, self.kp_signal))
            if do_it:
                msg  = "Model:\n sep=%6.2f +/- %5.2f mas" % \
                    (p[0], self.covar[0,0]*self.rms**2)
                msg += "\n  PA=%6.2f +/- %5.2f deg" % \
                    (p[1], self.covar[1,1]*self.rms**2)
                msg += "\n    c=%6.2f +/- %5.2f" % \
                    (p[2], self.covar[2,2]*self.rms**2)
            else:
                msg  = "Model:\n sep=%6.2f" % (p[0])
                msg += "\n   PA=%6.2f" % (p[1])
                msg += "\n  c=%6.2f" % (p[2])

            msg += "\n(rms= %.2f deg)" % (self.rms,)

            plt.text(0.0*mm, -0.75*mm, msg, 
                      bbox=dict(facecolor='white'), fontsize=14)

            msg = "Target: %s\nTelescope: %s\nWavelength = %.2f um" % (
                self.info['fname'], self.info['tel'], self.info['filter']*1e6)

            plt.text(-0.75*mm, 0.5*mm, msg,
                      bbox=dict(facecolor='white'), fontsize=14)
        except:
            pass

        plt.xlabel('Data Kernel-Phase Signal (deg)')
        plt.ylabel('Signal - Model Residuals (deg)')
        plt.draw()

    # =========================================================================
    # =========================================================================

    def extract_from_fits_frame(self, fname):
        ''' Extract the Kernel-phase signal from a single fits frame.

        ----------------------------------------------------------------
        Assumes that the fits file has been somewhat massaged, and is
        pretty much ready to go: no pairwise subtraction necessary, 
        etc...
        
        In addition to the actual Kernel-phase signal, some information
        is extracted from the fits header to help with the interpretation
        of the data to follow (think orientation of the telescope!).
        ----------------------------------------------------------------
        '''

        hdr = pf.getheader(fname)
        if hdr['TELESCOP'] == 'Keck II': self.info = get_keck_keywords(hdr)
        if hdr['TELESCOP'] == 'HST':     self.info = get_nic1_keywords(hdr)
        
        # read and fine-center the frame
        im = recenter(pf.getdata(fname), [80, 40, 10, 5], 25.0)

        sz = im.shape[0] # image is now square
        dz = sz/2.

        # meter to pixel conversion factor
        m2pix = mas2rad(self.info['pscale']) * sz / self.info['filter']
        uv_samp = self.uv * m2pix + dz # uv sample coordinates in pixels

        # calculate and normalize Fourier Transform
        ac = shift(fft(shift(im)))
        ac /= (np.abs(ac)).max() / self.nbh

        self.ac = ac
        self.im = im

        #for j in range(self.nbuv):
        xx = np.cast['int'](np.round(uv_samp[:,0]))
        yy = np.cast['int'](np.round(uv_samp[:,1]))
        data_cplx = ac[xx, yy]

        # ---------------------------
        # calculate the Kernel-phases
        # ---------------------------
        self.kp_signal = np.dot(self.KerPhi, 
                                self.RED*np.angle(data_cplx)) / dtor
        
        return self.kp_signal
        

    # =========================================================================
    # =========================================================================

    def plot_im_and_ac(self, rad=32):
        '''Plot of the image and its Fourier Transform.
        
        -----------------------------------------------------------------
        Include a plot of the sampled uv points overlaid with the data.
        Not always necessary, but a useful plot to look at sometimes...
        ----------------------------------------------------------------- '''

        sz = self.im.shape[0]
        dz = sz/2

        # meter to pixel conversion factor
        m2pix = mas2rad(self.info['pscale']) * sz / self.info['filter']
        uv_samp = self.uv * m2pix + dz # uv sample coordinates in pixels
        # pixel to mas conversion factor
        #pix2mas = self.info['pscale']/sz
        uvr = np.round(np.max(np.abs(uv_samp)), -1)/2

        plt.clf()
        f0 = plt.subplot(121)
        f0.imshow(self.im[dz-rad+1:dz+rad,dz-rad+1:dz+rad]**0.3)
        plt.title('Image')
        f1 = plt.subplot(122)
        f1.imshow((np.angle(self.ac)))
        plt.title('Fourier-phase')
        f1.plot(uv_samp[:,0], uv_samp[:,1], 'bo')
        f1.axis([dz-uvr,dz+uvr,dz-uvr,dz+uvr])
        plt.draw()
        plt.show()

    # =========================================================================
    # =========================================================================

    def save_image(self,fname,rad=32,fmt='png'):
        '''Saves the image as per plot_im'''

        sz = self.im.shape[0]
        dz = sz/2

        plt.clf()
        plt.imshow(self.im[dz-rad+1:dz+rad,dz-rad+1:dz+rad]**0.3)
        plt.title('Image')
        plt.savefig(fname+'.'+fmt,format=fmt)

    # =========================================================================
    # =========================================================================

    def save_to_file(self, file):
        ''' Export the KerPhaseData structure into a pickle file
        
        ----------------------------------------------------------------
        In addition to the KernelPhase data structure, 
        actual Kernel-phase signal is included, along with information
        related to the data: wavelength, instrument, 
        ---------------------------------------------------------------- '''
        try: data = {
            # -------------------
            # from KernelPhase
            # -------------------
            'name'   : self.name,
            'mask'   : self.mask,
            'uv'     : self.uv,
            'TFM'    : self.TFM,
            'KerPhi' : self.KerPhi,
            'RED'    : self.RED,
            # -------------------
            # additional information
            # -------------------
            'info'     : self.info,
            'comments' : self.comments,
            # ------------------------
            # actual Kernel-phase data
            # ------------------------
            'kp_signal' : self.kp_signal,
            'kp_error'  : self.kp_error}

        except:
            print("KerPhase_Relation data structure is incomplete")
            print("File %s wasn't saved!" % (file,))
            return None
        # -------------
        try: myf = open(file, "w")
        except:
            print("File %s cannot be created."+
                  " KerPhase_Relation data structure wasn't saved." % (file,))
            return None
        # -------------
        pickle.dump(data, myf)
        myf.close()
        
    # =========================================================================
    # =========================================================================

    def load_cals(self,calfname):


        '''Load the calibration specific information.

        If the file happens to be incomplete, the object will contain 
        a bunch of empty arrays.'''
    
        try:
            myf = open(calfname, "r")
            cals = pickle.load(myf)
            myf.close()
        except:
                print("File %s isn't a valid calibration structure" % (file))
                return None

        try:    self.kmean10879H = cals['kmean10879H']
        except:
            self.kmean10879H = 0.0
            print 'H band mean calibration failed to load'

        try:    self.kmean10879J = cals['kmean10879J']
        except:
            self.kmean10879J = 0.0
            print 'J band mean calibration failed to load'

        try:    self.krms10879J = cals['krms10879J']
        except:
            self.krms10879J = 0.0
            print 'J band RMS calibration failed to load'

        try:    self.krms10879H = cals['krms10879H']
        except:
            self.krms10879H = 0.0
            print 'H band RMS calibration failed to load'

        try:    self.kmean10143H = cals['kmean10143H']
        except:
            self.kmean10143H = 0.0
            print 'H band mean calibration failed to load'

        try:    self.kmean10143J = cals['kmean10143J']
        except:
            self.kmean10143J = 0.0
            print 'J band mean calibration failed to load'

        try:    self.krms10143J = cals['krms10143J']
        except:
            self.krms10879J = 0.0
            print 'J band RMS calibration failed to load'

        try:    self.krms10143H = cals['krms10143H']
        except:
            self.krms10143H = 0.0
            print 'H band RMS calibration failed to load'

    # =========================================================================
    # =========================================================================


    def calibrate(self,kmean,krms):

        '''Subtract kernel phase ensemble mean from signal and update
        signal and attribute properties accordingly. '''

##        if np.shape(kmean) is not ():
##            self.kp_signal = np.subtract(self.kp_signal,kmean)
##        elif np.shape(kmean) is():
##            self.kp_error = np.subtract(self.kp_signal,kmean*np.ones(self.KerPhi.shape[0]))
##        else: print "Something went horribly wrong in the systematic!"
##        
##        if np.shape(krms) is not ():
##            self.kp_signal = krms
##        elif np.shape(krms) is (): #scalar error
##            self.kp_error = krms*np.ones(self.KerPhi.shape[0])
##        else: print "Something went horribly wrong in assigning errors!!"

        self.kp_error = krms

        return None

    # =========================================================================
    # =========================================================================

    def mh(self,init_params,ranges,nsteps,sepmin=30.,sepmax=100.,
           anglemin=0.,anglemax=360.,cmin=1.01,cmax=10.):

        ''' A Metropolis-Hastings implementation for fitting kernel-phase data.
        Requires 3-element numpy arrays init_params and ranges, and a number of
        steps nsteps. Uses a Jeffreys prior for separation and contrast and a
        uniform prior for angle. This has an optional range determined by keywords'''

        params = init_params #initialise the parameter vector

        nsteps = int(nsteps)

        test = phase_binary(self.uv[:,0], self.uv[:,1],
                            self.info['filter'], [params[0],params[1]+90.0-self.info['orient'],params[2]])
        modl_ker = np.dot(self.KerPhi, self.RED*test)
        chi2 = np.sum(((modl_ker-self.kp_signal)/self.kp_error)**2)

        # work out prior odds factors

        prs = 1./np.log(sepmax/sepmin)
        pra = 1./(anglemax-anglemin)
        prc = 1./np.log(cmax/cmin)

        prior = 1/params[0]/params[2]*prs*pra*prc

        samples = np.zeros((3,nsteps))
        naccepted = 0
        nrejected = 0
        
        for j in range(0,int(nsteps)):

            #obtain a sample from our trial distribution

            new_params = params + np.random.randn(3)*ranges

            new_params[1] = np.mod(new_params[1],360.) # to avoid wrapping

            print 'Iteration',j, 'Parameters',new_params

            #compute the new chi2
            
            test = phase_binary(self.uv[:,0], self.uv[:,1],
                                self.info['filter'], [new_params[0],new_params[1]+90.0-self.info['orient'],new_params[2]])
            modl_ker = np.dot(self.KerPhi, self.RED*test)
            newchi2 = np.sum(((modl_ker-self.kp_signal)/self.kp_error)**2)
            print 'Chi-squared', newchi2

            # calculate the prior odds - Jeffreys for sep, c, uniform for angle

            newprior = 1/new_params[0]/new_params[2]*prs*pra*prc

            # work out the Metropolis ratio

            if sepmin<new_params[0]<sepmax and anglemin<new_params[1]<anglemax and cmin<new_params[2]<cmax:
                
                ratio = newprior/prior * np.exp(0.5*(chi2-newchi2))

            else: ratio = 0 #reject an out of range proposal

            # do we accept?

            u = np.random.rand()

            if ratio >= u:
                params = new_params
                prior = newprior
                chi2 = newchi2
                print 'Accepted',ratio
                naccepted += 1
            else:
                print 'Rejected',ratio
                nrejected +=1
            
            samples[:,j]= params

            #print 'Accepted',naccepted,'Rejected',nrejected
                
            print "Acceptance Fraction",naccepted/float(nrejected+naccepted)
        return samples
            
