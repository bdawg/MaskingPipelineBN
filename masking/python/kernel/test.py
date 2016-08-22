import pyker
import numpy as np
import pylab as plt
import pickle
from pyker_tools import *

x = np.linspace(0,10,100)
y = np.sin(x)

plt.plot(x,y)

plt.draw()
