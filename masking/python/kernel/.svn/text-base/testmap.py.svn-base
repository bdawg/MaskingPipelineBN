import pyker
from pyker_tools import *
import numpy as np, matplotlib.pyplot as plt
#a.fit_map(sep=524,cmin = 1.01,cmax=2.0)

ddir = './data/GO10143/'

# extract a

a = pyker.KernelPhaseData("./kerphi/hst.pick")
a.name = "HST - 10143 trial"
a.extract_from_fits_frame(ddir + 'n8yj16010_mos.fits')
a.load_cals('calibs.pick')
a.calibrate(a.kmean10143H,a.krms10143H)

#set up grid

nx, ny = 100,120
sep = 524
x = np.linspace(0,359.9,nx)
y = np.linspace(2.5,3.5,ny)
xx,yy  = np.meshgrid(x,y)

#calculate models

models = [[phase_binary(a.uv[:,0],a.uv[:,1],1.71e-6,
                        [524,xx[j][i]+90.0-a.info['orient'],yy[j][i]])for i in range(0,nx)] for j in range(0,ny)]
#calculate kerphases
kphasemodels = [[np.dot(a.KerPhi,a.RED*models[j][i]) for i in range (0,nx)]
                for j in range(0,ny)]

#calculate departures

diff = [[np.divide(np.subtract(a.kp_signal,kphasemodels[j][i])**2,
                   a.kp_error) for i in range(0,nx)] for j in range(0,ny)]

chi2 = [[np.sum(diff[j][i]) for i in range(0,nx)] for j in range(0,ny)]

chi2 = np.array(chi2)/(a.nkphi-4)

xlbl = 'Position Angle (degrees)'
ylbl = 'Contrast'

z = xx**2 + yy**2
clevs = chi2.min()*np.linspace(0,10,5)
plt.contourf(x,y,chi2,levels = clevs,cmap=plt.cm.bone)
plt.colorbar()
plt.xlabel(xlbl)
plt.ylabel(ylbl)
plt.show()
