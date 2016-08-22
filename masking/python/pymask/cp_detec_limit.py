import numpy as np
import matplotlib.pyplot as plt
import pickle
import time
import pymask
import multiprocessing

'''--------------------------------------------------------------
cp_detec_limit.py - uses a Monte Carlo simulation to establish
contrast-separation detection limits given an array of standard
deviations per kernel phase.

Written by F. Martinache and B. Pope.

Because different separation, contrast grid points are entirely
separate, this task is embarrassingly parallel. If you want to 
speed up the calculation, use multiprocessing with a threads 
argument equal to the number of available cores.

--------------------------------------------------------------'''
#ddir = '/home/ben/Research/Tim/'#
ddir = '/home/ben/Research/Desert/'

#oifitsname = 'MIRC_L2.Theta_Cyg.2012Jun19_JDM_2012Jun29.XCHAN.JDM.AVG15m'#
oifitsname = 'Binary1'#
#oifitsname = 'ScoCen2_Targ1'#

oifitsfile = oifitsname+'.oifits' 

addederror = 0.0

#------------------------
# first, load your data!
#------------------------

a = pymask.cpo(ddir+oifitsfile)

limits = pymask.detec_limits(a,threads=8,cmax=1000,nsep=64,ncon=64)
