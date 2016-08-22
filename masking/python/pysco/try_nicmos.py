#!/usr/bin/env python

import pysco
import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

# provide an absolute path to the data directory
# ----------------------------------------------
home = os.environ['HOME']
ddir = './data/' # data directory

''' ------------------------------------------------------------------
    For this demonstration, two HST data sets on the same target at two
    different wavelengths (110W and 170M) are provided as part of this
    initial package.
    The target is an obvious binary, so the final correlation plot 
    should be a satisfying one (especially at 1.7 um).
    
    after the "import pysco" command, the documentation for the code is
    available typing "help(pysco)", in the python command line.
    
    Cheers,
    
    Frantz.
    ----------------------------------------------------------------- '''

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
a.extract_KPD(ddir+'n8yj59020_mos.fits.gz', plotim=True)
a.kpi.name = "2M XXXX-XX" #  # labels the data

# ---------------------------------------
# optional plot you may want to look at
# ---------------------------------------
#plt.figure(2, (10,5))
#a.kpi.plot_pupil_and_uv(1.5) # optional plot of Kernel-phase model


# ------------------------
# 2. model fit of the data
# ------------------------
import pysco.fitting as fit

# binary_KPD_fit uses scipy's leastsq procedure, which is itself an 
# implementation of Levenberg-Marquardt Algorithm (LMA) to minimize the
# variance between the data and a binary-star model. Just like every
# implementation of LMA, you need to provide a initial set of
# parameters.
# it returns the best fit parameters, and a covariance matrix that can
# be used to determine uncertainties on the fit.

params0 = [150.0, 80.0, 5.0]    # initial parameters for model-fit
optim   = fit.binary_KPD_fit(a, params0)
params  = optim[0]    # best fit parameters (after least square)

# -------------------
# 3. correlation plot
# -------------------
fit.correlation_plot(a, params, plot_error=False)


#if optim[1] != None:
#    print "\nParmater estimate covariance matrix:"
#    print "-----------------------------------\n"
#    print np.round((a.rms)**2 * optim[1], 4)
