'''-----------------------------------------------------------------

kpstest.py - a script for using the output of kpsignal.py to mock up
a noisy phase-encoded signal and see what we can get out of it.

-----------------------------------------------------------------'''

import kpsignal
import numpy as np
import matplotlib.pyplot as plt

# set up the basic parameters

amplitude = 1.0 # carrier amplitude

n = 32 # number of channels - integer 32 is good

kp = kpsignal.KernelPhase(2*n+1) #initialise a 2n-bit (n channel) encoding

nkp = kp.nKPhi # create local variables for class attributes
TFM = kp.TFM
RED = kp.RED
KerPhi = kp.KerPhi
RowPhi = kp.RowPhi
KerPhiCol = kp.KPhiCol
RowPhiCol = kp.RPhiCol
nKPhi = kp.nKPhi
nRPhi = kp.nRPhi
V = kp.V
Vinv = kp.Vinv


# generate a random bit packet for the kernel space, encode as pi/0

Kbits = np.random.random_integers(0,1,size=n)*np.pi

# generate a random bit packet for the row space, encode as pi/0

Rbits = np.random.random_integers(0,1,size=n)*np.pi

# generate the singular phase projection of this bitstring

Sbits = np.zeros(2*n)

for j in range(0,nKPhi):
    Sbits[KerPhiCol[j]] = Kbits[j]

for j in range(0,nRPhi):
    Sbits[RowPhiCol[j]] = Rbits[j]

Sbits = np.matrix(Sbits)
Sbits = Sbits.T

# generate the frequency phase representation of this singular phase string

Frep = Vinv*Sbits

Frep = Frep.T[0]

Frep = [Frep[0,j] for j in range(0,Frep.shape[1])] # clumsy workaround to leave matrix notation

# now we generate a frequency spectrum

freqs = np.fft.fftfreq(2*n+1, d = 1e-6)

sideband = [np.complex(0,F) for F in Frep]

sideband = [np.exp(s) for s in sideband]

# now we generate the amplitudes

A = np.fft.irfft(sideband,axis = freqs)

Are = np.real(A)
Aim = np.imag(A)

plt.plot(Are**2 + Aim**2)
plt.show()
