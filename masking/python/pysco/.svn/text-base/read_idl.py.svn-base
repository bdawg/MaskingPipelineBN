import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle
import os
import sys
import gzip
import pysco
import re
import pysco.fitting as fit

from scipy.io.idl import readsav

ddir = '/home/ben/Research/Desert/'
mfdir = '/home/ben/code/masking/templates/'

'''-----------------------------------------------------------------------------
read_idl.py - reads in idlvar files for pysco
These idlvar files are read in as dictionaries with keys:

phs_v2corr
v2
cp
bs_var
v2_cov
ps
v2_cor
cp_cov
avar
bs_v2_cov
bs_cov
cp_var
err_avar
v2_all
bs
v2_sig
cvis_all
bs_all
cp_sig
mf_file - location of the matched filter file (in masking/templates/ [mf_file] )
u
v
bs_all: the bispectra from individual frames
-----------------------------------------------------------------------------'''

num = 460

filename = ddir+'bs0'+str(num)+'.idlvar'

data = readsav(filename)

#---------------------------------------
# now find a kpi appropriate to the mask
#---------------------------------------

# find appropriate mf files

mfdata = readsav(mfdir+data.mf_file)
filter = mfdata.filter[0]
maskfile = 'nirc2/info_g'+str(mfdata.n_holes)+'.idlvar'

#---------------------------------------
# Generate and save a KPI structure
#---------------------------------------

a = pysco.KPI(mfdir+data.mf_file,mfdir+maskfile)
a.name = 'g'+str(mfdata.n_holes)
a.save_to_file('./'+a.name+'.kpi.gz')

# -------------------
# 2. load the dataset
# -------------------

a = pysco.KPO('g9.kpi.gz')
a.extract_idl(filename,mfdir)
a.kpi.name = "2M XXXX-XX" #  # labels the data
a.filter = filter
a.vis_error = 0.05*np.ones(a.kpi.nbuv)

# ------------------------
# 3. model fit of the data
# ------------------------

npoints = 100 # active points
nsteps = 2400 # iterations
njumps = 20 # random walk steps

sepmin,sepmax,anglemin,anglemax,cminh,cmaxh,cminj,cmaxj = 100.,160.,0.,360.,1.00001,4.0,1.00001,4.0
cmin = 1.0001
cmax = 4.0

mode = 'kp' #set to 'kp' for kernel phase and 'vis' for vis2 fitting


[sepmean,anglemean,cmeanh,cmeanj,dsep,dangle,dch,dcj,logE,dlogE] = fit.nested(a,prange=[sepmin,sepmax,anglemin,anglemax,cmin,cmax],
                                                                 nsteps=nsteps,npoints=npoints,njumps=njumps,mode=mode)
