#!/usr/bin/env python

import pysco
import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
import pysco.fitting as fit

# provide an absolute path to the data directory
# ----------------------------------------------
home = os.environ['HOME']
ddir = './data/' # data directory


''' ------------------------------------------------------------------
    This is a demonstration of single-band nested sampler data analysis
    with pysco. Enter the appropriate prior range below, and choose the
    algorithm parameters, and it should work just fine.

    Note that vis2 fitting isn't quite up to speed yet, couple of bugs to
    iron out still. 

    Cheers,

    Ben
    ----------------------------------------------------------------- '''

npoints = 100 # active points
nsteps = 2400 # iterations
njumps = 20 # random walk steps

sepmin,sepmax,anglemin,anglemax,cminh,cmaxh,cminj,cmaxj = 100.,160.,0.,360.,1.00001,4.0,1.00001,4.0
cmin = 1.0001
cmax = 4.0

mode = 'kp' #set to 'kp' for kernel phase and 'vis' for vis2 fitting

# -------------------------------
# 1. create the KP info structure
# -------------------------------

# once saved, the kpi.gz structure can be directly reloaded when 
# creating a KPO instance, such as done in step #2.

a = pysco.KPI("./geometry/nicmos.txt")
a.name = "HST - NIC1" #  # add a label to the template
a.save_to_file('./hst.kpi.gz')

# -------------------
# 2. load the dataset
# -------------------

# load the FITS frame, and extract the Kernel-phases using the
# HST KPI template calculated at the previous step. 
# Two data sets are provided:
# n8yj59010_mos.fits.gz and 'n8yj59020_mos.fits.gz

a = pysco.KPO('./hst.kpi.gz')
a.extract_KPD(ddir+'n8yj59020_mos.fits', plotim=True)
a.kpi.name = "2M XXXX-XX" #  # labels the data
a.vis_error = 0.05*np.ones(a.kpi.nbuv)

c = pysco.KPO('./hst.kpi.gz')
c.extract_KPD(ddir+'n8yj59010_mos.fits', plotim=True)
c.kpi.name = "2M XXXX-XX" #  # labels the data
c.vis_error = 0.05*np.ones(a.kpi.nbuv)
d = [a,c]

# ---------------------------------------
# optional plot you may want to look at
# ---------------------------------------
#plt.figure(2, (10,5))
#a.kpi.plot_pupil_and_uv(1.5) # optional plot of Kernel-phase model


# ------------------------
# 3. model fit of the data
# ------------------------

[sepmean,anglemean,cmeanh,cmeanj,dsep,dangle,dch,dcj,logE,dlogE] = fit.nested(d,prange=[sepmin,sepmax,anglemin,anglemax,cminh,cmaxh,cminj,cmaxj],
                                                                 nsteps=nsteps,npoints=npoints,njumps=njumps,mode=mode)

# ------------------------
# 4. make pretty pictures
# ------------------------

params = [sepmean,anglemean,cmeanh,cmeanj]

b = fit.binary_KPD_model(a,params)
plt.figure(3)
fit.correlation_plot(a,b=b)

plt.show()
