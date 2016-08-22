#!/usr/bin/env python

import pyker
import numpy as np
import matplotlib.pyplot as plt

ddir = './data/' # data directory

''' ------------------------------------------------------------------
For this demonstration, two HST data sets on the same target at two
different wavelengths (110W and 170M) are provided as part of this
initial package. 
The target is an obvious binary, so the final correlation plot should
be a satisfying one (especially at 1.7 um).

after the "import pyker" command, the documentation for the code is
available typing "help(pyker)".

----------------------------------------------------------------- '''


# -------------------
# 1. load the dataset
# -------------------

# create an instance of KernelPhaseData using the HST template.
a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.load_cals("calibs.fits")
print "about to name"
a.name = "HST - NIC1" #  # add a label to the template

# load the FITS frame, and extract the Kernel-phases using the already
# loaded HST template. Two data sets are provided:
# n8yj59010_mos.fits.gz and 'n8yj59020_mos.fits.gz
print "about to extract"
a.extract_from_fits_frame(ddir + 'n8yj59020_mos.fits.gz')

# --------------------------------------------------------
# a couple of optional plots which you may want to look at
# --------------------------------------------------------
#a.plot_pupil_and_uv(1.5) # optional plot of Kernel-phase model
a.plot_im_and_ac()       # optional plot of image and Fourier counterpart

# ------------------------
# 2. model fit of the data
# ------------------------

# binary_fit uses Levenberg-Marquardt Algorithm (LMA) to minimize the
# variance between the data and a binary-star model. Just like every
# implementation of LMA, you need to provide a initial set of
# parameters.
# it returns the best fit parameters, and a covariance matrix that can
# be used to determine uncertainties on the fit.
print "up to model fitting"
init_params = [50.0, 40.0, 1.0]    # initial parameters for model-fit
optim = a.binary_fit(init_params)  # fit data with binary model
params = optim[0]                  # best fit parameters

b = pyker.KernelPhaseData("./kerphi/hst.pick")
b.binary_model(params, a.info)

# -------------------
# 3. correlation plot
# -------------------
a.correlation_plot(b)

if optim[1] != None:
    print "\nParameter estimate covariance matrix:"
    print "-----------------------------------\n"
    print np.round((a.rms)**2 * optim[1], 4)

if optim[0] != None:
    print "\nParameter estimate:"
    print "-----------------------------------\n"
    print optim[0]
plt.show()
