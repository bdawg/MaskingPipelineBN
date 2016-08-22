#!/usr/bin/python


#;This function takes a wavelength vector (in microns) and a flux
#;vector in Jy and returns a vector of magnitudes. It linearly
#;interpolates between the bands, and uses values from Allen's
#;Astrophysicsl quantities
#;The option of filter = 'J'|'H'|'K'|'V'|'R'|'I' just finds the mean
#;flux in the filter bandwidth, as defined roughly by Scholz model files
#;and Allen.

#
# converted to Python
# PNS Jun 2012
#
#
#
#
# USAGE
#   from jymag_conv import *
# then run as either
#   jy2mag(wavelengths, flux)
#   mag2jy(wavelengths, mag)
#


def jy2mag( wavelengths, flux, filt=''):
    """jy2mag( wavelengths, flux, filt='')"""
    from scipy.interpolate import interp1d
    from numpy import log10, array, ones
    wavelengths = array(wavelengths)
    flux = array(flux)
    #;Values for a 0 mag star. Note that I've doctored the B filter a
    #;little for cubic spline interpolation, maitaining the mean flux
    #;within the filter.
    filt_l =  array([0.36,0.41,0.47,0.55, 0.7,0.9, 1.22,1.63,2.19,3.45,  5,  10.2,21.0])
    filt_jy = array([1840,4240,4320,3832,2842,2246,1570,1020,636,313,182.5,42.7,10.0])#;J was 1724
    zmag_flux = (interp1d(filt_l,filt_jy, kind="cubic"))(wavelengths)
    
    if(wavelengths.max() > filt_l.max()) or (wavelengths.min() < filt_l.min()):
        print '### ERROR - Lambda out of interpolation range ### jy2mag(lambdas,fluxs)'
    
    ################
    #this part may not have survived the python conversion intact...
    # needs further testing
    if filt != '' and wavelengths.size != 1 :
        w = [-1]
        if (filt == 'J'):           w = (wavelengths >= 1.1  ) * (wavelengths <= 1.4)
        if (filt == 'H'):           w = (wavelengths >= 1.483) * (wavelengths <= 1.783)
        if (filt == 'K'):           w = (wavelengths >= 1.995) * (wavelengths <= 2.395)
        if (filt == 'L'):           w = (wavelengths >= 3.04 ) * (wavelengths <= 3.84)
        if (filt == 'L*'):          w = (wavelengths >= 3.5  ) * (wavelengths <= 4.1)
        if (filt == 'V'):           w = (wavelengths >= 0.505) * (wavelengths <= 0.595)
        if (filt == 'R'):           w = (wavelengths >= 0.595) * (wavelengths <= 0.805)
        if (filt == 'I'):           w = (wavelengths >= 0.8  ) * (wavelengths <= 1.0)
        if (filt == 'B'):           w = (wavelengths >= 0.39 ) * (wavelengths <= 0.49)
        if (filt == '78'):          w = (wavelengths >= 0.7773)* (wavelengths <= 0.7863)
        if (filt == '87'):          w = (wavelengths >= 0.8736)* (wavelengths <= 0.8818)
        if (filt == '70'):          w = (wavelengths >= 0.695) * (wavelengths <= 0.705)
        if (filt == '71'):          w = (wavelengths >= 0.710) * (wavelengths <= 0.720)
        if (filt == '905 box'):     w = (wavelengths >= 0.878) * (wavelengths <= 0.930)
        if (filt == '1040 narrow'): w = (wavelengths >= 1.035) * (wavelengths <= 1.045)
	 
        nw = w.size
        if nw != 1:
            mns = 0.5*(flux[w[0:nw-2]] + flux[w[1:nw-1]])
            diffs = wavelengths[w[1:nw-1]] - wavelengths[w[0:nw-2]]
            mnl = 0.5*(wavelengths[w[1:nw-1]] +wavelengths[w[0:nw-2]])
            weights = ones(mns.size) #;an array of the same size...
            #if (filter == 'L') then weights = interpol([0,0.5,0.85,0.86,0.19,0],[3.04,3.24,3.48,3.68,3.76,3.84],mnl)
            #if (filter == 'V') then weights = interpol([0,0.78,1.0,0.898,0.359,0.08,0],[0.48,0.51,0.53,0.55,0.6,0.64,0.68],mnl)
            mnflux = (mns*weights*diffs/mnl**2).sum()/(diffs/mnl**2).sum()
            #;Now do the same for zmag_flux
            zmag_flux = 0.5*(zmag_flux[w[0:nw-2]] + zmag_flux[w[1:nw-1]])
            mnzflux = (zmag_flux*weights*diffs/mnl**2).sum()/(diffs/mnl**2).sum()
            return 2.5*log10(mnzflux/mnflux)
        else:
            print 'Error: No input wavelengths within filter...'
            return 0
	#
    # The above needs more testing.
    ##############

    return 2.5*log10(zmag_flux/flux)

def mag2jy(wavelengths,mag_vector):
	"""mag2jy(wavelengths,mag_vector)"""
	from scipy.interpolate import interp1d
	from numpy import array
	wavelengths = array(wavelengths)
	mag_vector = array(mag_vector)

	#;bands=  U, B, V, R, I, J, H, K, L, M, N
	#lambd=[0.36, 0.44, 0.55, 0.70, 0.90, 1.25, 1.65, 2.2, 3.4, 5.0, 10.2]
	#jansky0=[1880., 4440, 3810, 2880, 2240, 1771, 1062, 629, 312, 180, 43]

	#; NEW set of standards from VEGA (Allen, astrophysical quantities)
	#; _EXCEPT_ ones indicated in brackets taken from above
	#; bands= (U),(B),V,(R),(I),J,H,Ks,K,L,L',M,8.7,N,11.7,Q

	lambd=  array([0.36,0.44,0.5556,0.70,0.90,1.215,1.654,2.157,2.179,3.547, 3.761,4.769,8.756,10.472,11.653,20.130])
	jansky0=array([1880,4440,3540,  2880,2240,1630, 1050, 667,  655,  276,  248,  160,  50.0, 35.2,  28.6,  9.70  ])

	if(wavelengths.max() > lambd.max()) or (wavelengths.min() < lambd.min()):
	    print '### ERROR - Lambda out of interpolation range ### mag2jy(lambdas,mags)'
        
	#jansky=spline(lambd,jansky0,wavelengths)
	jansky = (interp1d(lambd,jansky0, kind="cubic"))(wavelengths)
	#print mag_vector/2.5
	#print jansky
	#; convert magnitudes
	star_flx  = 10**(mag_vector/(-2.5)) * jansky

	return star_flx
