#plot kernel phase histograms

import pyker
import numpy as np
import matplotlib.pyplot as plt

plt.clf()

f1 = plt.hist(kmean,30)
f1.xlabel('Kernel Phase Signal Ensemble Mean (deg)')
f1.ylabel('Occurrences')
f1.title('Histogram of kernel phase mean values')
f2 = plt.hist(krms,30)
f2.xlabel('Kernel Phase Ensemble RMS')
f2.ylabel('Occurrences')
f2.title('Histogram of kernel phase mean values')
plt.draw()
plt.show()
