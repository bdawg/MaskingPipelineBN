'''---------------------------------------------------------------

kpsignal.py - a piece of code designed to test whether the kernel-
phase basis is useful for signal processing. We generate a random
bit string and encode this as 0 or pi kernel phases and then assign
row phases for convenience.

We then generate a phase vector from these by inverting the transfer
and apply these phases to a uniform-intensity carrier wave.

To receive the signal you interfere this with a local oscillator to
receive the GHz-freqency fringes. Then you FFT the whole fringe packet
and from its frequency representation you reconstruct the kernel and
row phases.

We will use a piece of code cannibalised from Frantz' KernelPhase
class object in pyker.py to generate the kernel and row phases from
a specified list of times.

---------------------------------------------------------------'''

import numpy as np, matplotlib.pyplot as plt, pyfits as pf
import pickle, os, pdb, matplotlib.cm as cm
from scipy.special import j1

from pyker_tools import *

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

dtor = np.pi/180.0

class KernelPhase(object):
    ''' Fundamental kernel-phase relations

    -----------------------------------------------------------------------
    This object condenses all the knowledge about a given instrument pupil 
    geometry into a series of arrays useful for kernel-phase analysis as 
    well as for other purposes, such as wavefront sensing.
    ----------------------------------------------------------------------- '''

    def __init__(self,n):
        try: self.from_coords(n)
        except: "You've done it wrong!"
    
    # =========================================================================
    # =========================================================================

    def from_coords(self,n):
        ''' Creation of the KerPhase_Relation object from a pupil mask file:

        ----------------------------------------------------------------
        This is the core function of this class, really...

        Input is a pupil coordinates file, containing one set of (x,y) 
        coordinates per line. Coordinates are in meters. From this, all 
        the intermediate products that lead to the kernel-phase matrix 
        KerPhi are calculated.
        ---------------------------------------------------------------- '''
        #self.mask = 1.0 * np.loadtxt(file) # sub-Ap. coordinate files
        self.mask = np.linspace(0,1,n) #linear time samples over 1 us
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

        k = 0 # index for possible combinations (k = f(i,j))
        
        uvi = np.zeros(nbh * (nbh-1), dtype=int) # arrays to store the possible
        uvj = np.zeros(nbh * (nbh-1), dtype=int) # combinations k=f(i,j) !!

        for i in range(nbh):     # do all the possible combinations of
            for j in range(nbh): # sub-apertures
                if i != j:
                    uvx[k] = self.mask[i] - self.mask[j]
                    uvi[k], uvj[k] = i, j
                    k+=1

        a = np.unique(np.round(uvx, ndgt)) # distinct u-component of baselines
        nbx    = a.shape[0]                # number of distinct u-components
        #uv_sel = np.zeros((0,2))           # array for "selected" baselines
        print nbx

        self.nbuv = np.shape(a)[0]/2 # actual number of distinct uv points
        self.uv   = a[:self.nbuv]  # discard second half (symmetric)
        print "%d distinct baselines were identified" % (self.nbuv,)

        # 2. Calculate the transfer matrix and the redundancy vector
        # --------------------------------------------------------------
        self.TFM = np.zeros((self.nbuv, self.nbh), dtype=float) # matrix
        self.RED = np.zeros(self.nbuv, dtype=float)             # Redundancy

        for i in range(self.nbuv):
            a=np.where(np.abs(self.uv[i]-uvx) <= prec)
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
        self.KPhiCol     = np.where(abs(S1) < 1e-3)[0]
        self.nRPhi  = np.size(np.where(abs(S1) > 1e-3)) # number of Row-phases
        self.RPhiCol     = np.where(abs(S1) > 1e-3)[0]
        self.KerPhi = np.zeros((self.nKPhi, self.nbuv)) # allocate the arrays
        self.RowPhi = np.zeros((self.nRPhi, self.nbuv))
        self.V      = Vh
        self.Vinv   = np.linalg.inv(Vh)

        for i in range(self.nKPhi):
            self.KerPhi[i,:] = (Vh)[self.KPhiCol[i],:]

        for i in range(self.nRPhi):
            self.RowPhi[i,:] = (Vh)[self.RPhiCol[i],:]

        #print '-------------------------------'
        #print 'Singular values for this array:\n', np.round(S, ndgt)
        #print '\nRedundancy Vector:\n', self.RED
        #self.name = array_name

    # =========================================================================
    # =========================================================================
