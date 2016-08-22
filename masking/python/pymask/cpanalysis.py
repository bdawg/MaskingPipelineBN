import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf
import pickle
import gzip
import pymultinest
import os, threading, subprocess
import matplotlib.pyplot as plt
import json
import time

import pymask

'''---------------------------------------------------------
cpanalysis.py - a script for analysing masking data using 
either the MCMC Hammer ensemble affine invariant MCMC 
algorithm, or the MultiNest multimoded nested sampling 
algorithm. 

Comment out whichever you don't want!

This depends significantly on emcee, the Python MCMC Hammer,
and PyMultiNest, a wrapper for the Fortran-based MultiNest.
---------------------------------------------------------'''

ddir = '/home/ben/Research/Desert/'

oifitsfile = 'Binary1.oifits' #'ScoCen2_Targ1.oifits'

#------------------------
# first, load your data!
#------------------------

a = pymask.cpo(ddir+oifitsfile)

#--------------------------
# Use this for MCMC Hammer
# Warning - get ivar right!
#--------------------------

# ivar = [52., 192., 1.53] # this has to be pretty much near the peak or MCMC Hammer can get lost

# pymask.hammer(a,ivar=ivar,plot=True)

#--------------------------
# Use this for MultiNest
#--------------------------

paramlimits = [20.,250.,0.,360.,1.0001,2] #sepmin,sepmax,anglemin,anglemax,cmin,cmax

model = pymask.nest(a,paramlimits=paramlimits,multi=False)