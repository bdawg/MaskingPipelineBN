''' --------------------------------------------------------------------
                PYSCO: PYthon Self Calibrating Observables
    --------------------------------------------------------------------
    ---
    pysco is a python module to create, and extract Kernel-phase data 
    structures, using the theory of Martinache, 2010, ApJ, 724, 464.
    ----

    This file contains the definition of the KPO class:
    --------------------------------------------------

    an object that contains Ker-phase information (kpi), data (kpd) 
    and relevant additional information extracted from the fits header
    (hdr)
    -------------------------------------------------------------------- '''

#!/usr/bin/env python

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

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

from scipy.interpolate import griddata

from scipy.io.idl import readsav

import oifits

import core
from core import *

import kpi
from kpi import *

class KPO():
    ''' Class used to manipulate multiple Ker-phase datasets

        -------------------------------------------------------------------
        The class is designed to handle a single or multiple frames
        that can be combined for statistics purpose into a single data
        set.
        ------------------------------------------------------------------- '''

    def __init__(self, kp_fname):
        # Default instantiation.
        self.kpi = KPI(kp_fname)

        # if the file is a complete (kpi + kpd) structure
        # additional data can be loaded.
        try:
            myf = gzip.GzipFile(kp_fname, "r")
            data = pickle.load(myf)
            myf.close()

            self.kpd = data['kpd']
            self.kpe = data['kpe']
            self.hdr = data['hdr']
        except:
            print("File %s contains KPI information only" % (kp_fname,))

    # =========================================================================
    # =========================================================================
    def extract_KPD(self, path, plotim=False, ave="none"):
        ''' extract kernel-phase data from one or more files (use regexp).

        If the path leads to a fits data cube, or to multiple single frame
        files, the extracted kernel-phases are consolidated into a
        unique kpd object.      
        
        '''
        fnames = glob.glob(path)
        nf = fnames.__len__()
        if nf != 0:
            fits_hdr = pf.getheader(fnames[0])
        elif nf == 0:
            fits_hdr = pf.getheader(fnames)
        else:
            print 'Frame number error'
        
        print "%d frames will be open" % (nf,)

        hdrs = [] # empty list of kp data headers

        # =========================
        if fits_hdr['NAXIS'] < 3:
            kpds = np.zeros((nf, self.kpi.nkphi)) # empty 2D array of kp
            vis2s = np.zeros((nf,self.kpi.nbuv))
        
            for i, fname in enumerate(fnames):
                (hdr, sgnl, vis2, im, ac) = \
                    extract_from_fits_frame(fname, self.kpi, save_im=True)
                kpds[i] = sgnl
                vis2s[i] = vis2
                hdrs.append(hdr)
            if nf == 1:
                hdrs = hdr
                kpds = sgnl
                vis2s = vis2
                
        # =========================
        if fits_hdr['NAXIS'] == 3:
            kpds = np.zeros((fits_hdr['NAXIS3'], 
                             self.kpi.nkphi)) # empty kp array
            dcube = pf.getdata(fnames[0])
            nslices = fits_hdr['NAXIS3']
            for i in xrange(nslices):
                sys.stdout.write(
                    "\rextracting kp from img %3d/%3d" % (i+1,nslices))
                sys.stdout.flush()
                (hdr, sgnl, vis2) = extract_from_array(dcube[i], fits_hdr, self.kpi, 
                                                 save_im=False, re_center=True,
                                                 wrad=50.0, plotim=True)
                kpds[i] = sgnl
                vis2s[i]= vis2
                hdrs.append(hdr)

        self.kpe = np.std(kpds, 0)
        self.vis2e = np.std(vis2s,0)

        if ave == "median":
            print "median average"
            self.kpd = np.median(kpds, 0)
            self.vis2 = np.median(vis2s,0)
            self.hdr = hdrs[0]

        if ave == "mean":
            print "mean average"
            self.kpd = np.mean(kpds, 0)
            self.vis2 = np.mean(vis2s,0)
            self.hdr = hdrs[0]

        if ave == "none":
            print "no average"
            self.kpd = kpds
            self.vis2 = vis2s
            self.hdr = hdrs

    # =========================================================================
    # =========================================================================
    def extract_idl(self, path, mfpath, plotim=False, ave="none"):
        '''Extracts kernel phase data from an idlvar file.        
        '''

        data = readsav(path)

        self.hdr = get_idl_keywords(mfpath+data.mf_file)
        self.bispectrum = data.bs

        self.kpd = data.cp
        self.vis2 = data.v2

        self.kpe = data.cp_sig
        self.vis2e = data.v2_sig

        # self.kpe = np.std(kpds, 0)
        # self.vis2e = np.std(vis2s,0)

        # if ave == "median":
        #     print "median average"
        #     self.kpd = np.median(kpds, 0)
        #     self.vis2 = np.median(vis2s,0)
        #     self.hdr = hdrs[0]

        # if ave == "mean":
        #     print "mean average"
        #     self.kpd = np.mean(self.kpds, 0)
        #     self.vis2 = np.mean(self.vis2s,0)
        #     self.hdr = hdrs[0]

        # if ave == "none":
        #     print "no average"
        #     self.kpd = self.kpds
        #     self.vis2 = self.vis2s
        #     self.hdr = self.hdrs

    # =========================================================================
    # =========================================================================
    def copy(self):
        ''' Returns a deep copy of the Multi_KPD object.
        '''
        res = copy.deepcopy(self)
        return res

    # =========================================================================
    # =========================================================================
    def calibrate(self, calib, regul="None"):
        ''' Returns a new instance of Multi_KPD object.

        Kernel-phases are calibrated by the calibrator passed as parameter.
        Assumes for now that the original object and the calibrator are
        collapsed into one single kp data set. '''

        res = copy.deepcopy(self)

        if np.size(calib.kpd.shape) == 1:
            res.kpd -= calib.kpd
            return res
        else:
            coeffs = super_cal_coeffs(self, calib, regul)
            
        return res

    # =========================================================================
    # =========================================================================
    def average_KPD(self, algo="median"):
        ''' Averages the multiple KP data into a single series.

        Default is "median". Other option is "mean".
        '''
        
        temp = np.array(self.kpd)

        if algo == "median":
            aver = np.median(self.kpd, 0)
        else:
            aver = np.mean(self.kpd, 0)


        # ----------------------
        # update data structures
        # ----------------------
        self.kpd = aver

        # -------------------------------------
        # update data header (mean orientation)
        # -------------------------------------
        nh = self.hdr.__len__()
        ori = np.zeros(nh)

        for i in xrange(nh):
            ori[i] = self.hdr[i]['orient']

        self.hdr = self.hdr[0] # only keep one header
        self.hdr['orient'] = np.mean(ori)

        return self.kpd

    # =========================================================================
    # =========================================================================
    def save_to_file(self, fname):
        '''Saves the kpi and kpd data structures in a pickle

        --------------------------------------------------------------
        The data can then later be reloaded for additional analysis,
        without having to go through the sometimes time-consuming
        extraction from the original fits files.

        To save on disk space, this procedure uses the gzip module.
        While there is no requirement for a specific extension for the
        file, I would recommend that one uses ".kpd.gz", so as to make
        it obvious that the file is a gzipped kpd data structure.
        --------------------------------------------------------------
        '''

        try:
            data = {'name'   : self.kpi.name,
                    'mask'   : self.kpi.mask,
                    'uv'     : self.kpi.uv,
                    'TFM'    : self.kpi.TFM,
                    'KerPhi' : self.kpi.KerPhi,
                    'RED'    : self.kpi.RED}
        except:
            print("KPI data structure is incomplete")
            print("File %s was not saved to disk" % (fname,))
            return(None)

        try:
            data['hdr'] = self.hdr
            data['kpd'] = self.kpd
            data['kpe'] = self.kpd

        except:
            print("KPD data structure is incomplete")
            print("File %s was nevertheless saved to disk" % (fname,))

        try:
            myf = gzip.GzipFile(fname, "wb")
        except:
            print("File %s cannot be created" % (fname,))
            print("Data was not saved to disk")
            return(None)

        pickle.dump(data, myf, -1)
        myf.close()
        print("File %s was successfully written to disk" % (fname,))
        return(0)

    # =========================================================================
    # =========================================================================
    def plot_uv_phase_map(self, data=None, reso=400):

        uv = self.kpi.uv
        
        Kinv = np.linalg.pinv(self.kpi.KerPhi)

        dxy = np.max(np.abs(uv))
        xi = np.linspace(-dxy, dxy, reso)
        yi = np.linspace(-dxy, dxy, reso)

        if data == None:
            data = np.dot(Kinv, self.kpd)
        z1 = griddata((np.array([uv[:,0], -uv[:,0]]).flatten(),
                       np.array([uv[:,1], -uv[:,1]]).flatten()),
                      np.array([data, -data]).flatten(),
                      (xi[None,:], yi[:,None]), method='linear')
        
        z2 = griddata((np.array([uv[:,0], -uv[:,0]]).flatten(), 
                       np.array([uv[:,1], -uv[:,1]]).flatten()), 
                      np.array([self.kpi.RED, self.kpi.RED]).flatten(),
                      (xi[None,:], yi[:,None]), method='linear')

        plt.imshow(z1)
        return (z1, z2)
