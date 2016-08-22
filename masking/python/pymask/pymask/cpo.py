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

from cp_tools import *

'''------------------------------------------------------------------------
cpo.py - Python class for manipulating oifits format closure phase data.
------------------------------------------------------------------------'''


class cpo():
    ''' Class used to manipulate multiple closure phase datasets'''

    def __init__(self, oifits):
        # Default instantiation.

        # if the file is a complete (kpi + kpd) structure
        # additional data can be loaded.
        try:
           self.extract_from_oifits(oifits)
        except:
            print('Invalid file.')


    def extract_from_oifits(self,filename):
        '''Extract closure phase data from an oifits file.'''

        data = oifits.open(filename)
        self.name = ''

        self.ndata = len(data.t3)

        for j in data.wavelength:
            wavel = data.wavelength[j].eff_wave
            self.wavel = wavel[0]
            break

        self.target = data.target[0].target

        t3data = []
        t3err = []
        self.u = np.zeros((self.ndata,3))
        self.v = np.zeros((self.ndata,3))

        for j, t3 in enumerate(data.t3):
            t3data.append(t3.t3phi[0])
            t3err.append(t3.t3phierr[0])
            self.u[j,:] = [t3.u1coord,t3.u2coord,-(t3.u1coord+t3.u2coord)]
            self.v[j,:] = [t3.v1coord,t3.v2coord,-(t3.v1coord+t3.v2coord)]

        self.t3data = np.array(t3data)
        self.t3err = np.array(t3err)
