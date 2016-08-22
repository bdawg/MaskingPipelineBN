import pysco
import pysco.fitting as fit
import matplotlib.pyplot as plt
import time
import os
import pdb

# provide an absolute path to the data directory
# ----------------------------------------------

home = os.environ['HOME']
ddir = home+"/data/P3K_AOph/" # data directory
a = pysco.KPI("./geometry/palomar_med_cross.txt")
a.save_to_file('./palomar.pick.gz')

plt.figure(1)

a = pysco.KPO("palomar.kpi.gz")
b = pysco.KPO("palomar.kpi.gz")

a.extract_KPD(ddir + "cubeF_2056.fits", plotim=True, ave="median")
#a.extract_KPD(ddir + "bin5_N.fits",     plotim=True, ave="median")
b.extract_KPD(ddir + "cubeF_2881.fits", plotim=True, ave="median")

# calibrated kernel-phase data
# ----------------------------
c = a.copy() 
c.kpd = a.kpd - b.kpd
c.kpe = np.sqrt(a.kpe**2 + b.kpe**2)

plt.figure(2)

map = fit.chi2map_sep_pa(c, con=30, reduced=True)

params0 = [131.0, 276.0, 30.0]    # initial parameters for model-fit
optim   = fit.binary_KPD_fit(c, params0)
params  = optim[0]    # best fit parameters (after least square)

fit.correlation_plot(c, params)#, plot_error=False)
