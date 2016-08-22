import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.signal import medfilt2d as medfilt
from scipy.io.idl import readsav
import pyfits as pf

''' ================================================================
    small tools and functions useful for the manipulation of
    Ker-phase data.
    ================================================================ '''

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

dtor = np.pi/180.0

# =========================================================================
# =========================================================================

def mas2rad(x):
    ''' Convenient little function to convert milliarcsec to radians '''
    return x*np.pi/(180*3600*1000)

# =========================================================================
# =========================================================================

def rad2mas(x):
    ''' Convenient little function to convert radians to milliarcseconds '''
    return x/np.pi*(180*3600*1000)
# =========================================================================
# =========================================================================

def phase_binary(u, v, wavel, p):
    ''' Calculate the phases observed by an array on a binary star
    ----------------------------------------------------------------
    p: 3-component vector (+2 optional), the binary "parameters":
    - p[0] = sep (mas)
    - p[1] = PA (deg) E of N.
    - p[2] = contrast ratio (primary/secondary)
    
    optional:
    - p[3] = angular size of primary (mas)
    - p[4] = angular size of secondary (mas)

    - u,v: baseline coordinates (meters)
    - wavel: wavelength (meters)
    ---------------------------------------------------------------- '''

    p = np.array(p)
    # relative locations
    th = (p[1] + 90.0) * np.pi / 180.0
    ddec =  mas2rad(p[0] * np.sin(th))
    dra  = -mas2rad(p[0] * np.cos(th))

    # baselines into number of wavelength
    x = np.sqrt(u*u+v*v)/wavel

    # decompose into two "luminosity"
    l2 = 1. / (p[2] + 1)
    l1 = 1 - l2
    
    # phase-factor
    phi = np.zeros(u.size, dtype=complex)
    phi.real = np.cos(-2*np.pi*(u*dra + v*ddec)/wavel)
    phi.imag = np.sin(-2*np.pi*(u*dra + v*ddec)/wavel)

    # optional effect of resolved individual sources
    if p.size == 5:
        th1, th2 = mas2rad(p[3]), mas2rad(p[4])
        v1 = 2*j1(np.pi*th1*x)/(np.pi*th1*x)
        v2 = 2*j1(np.pi*th2*x)/(np.pi*th2*x)
    else:
        v1 = np.ones(u.size)
        v2 = np.ones(u.size)

    cvis = l1 * v1 + l2 * v2 * phi
    phase = np.angle(cvis, deg=True)
    return np.mod(phase + 10980., 360.) - 180.0

# =========================================================================
# =========================================================================

def vis2_binary(u, v, wavel, p):
    ''' -------------------------------------------------------------------
      Calculate the vis-squareds observed by an array on a binary star
    p: 3-component vector (+2 optional), the binary "parameters":
      p[0] = sep (mas)
      p[1] = PA (deg) E of N.
      p[2] = contrast ratio (primary/secondary)

    optional:
      p[3] = angular size of primary (mas)
      p[4] = angular size of secondary (mas)

    u,v: baseline coordinates (meters)
    wavel: wavelength (meters)
    ---------------------------------------------------------------- '''
    p = np.array(p)
    # relative locations
    th = (p[1] + 90.0) * np.pi / 180.0
    ddec =  mas2rad(p[0] * np.sin(th))
    dra  = -mas2rad(p[0] * np.cos(th))

    # baselines into number of wavelength
    x = np.sqrt(u*u+v*v)/wavel

    # decompose into two "luminosity"
    l2 = 1. / (p[2] + 1)
    l1 = 1 - l2
    
    # phase-factor
    phi = np.zeros(u.size, dtype=complex)
    phi.real = np.cos(-2*np.pi*(u*dra + v*ddec)/wavel)
    phi.imag = np.sin(-2*np.pi*(u*dra + v*ddec)/wavel)

    # optional effect of resolved individual sources
    if p.size == 5:
        th1, th2 = mas2rad(p[3]), mas2rad(p[4])
        v1 = 2*j1(np.pi*th1*x)/(np.pi*th1*x)
        v2 = 2*j1(np.pi*th2*x)/(np.pi*th2*x)
    else:
        v1 = np.ones(u.size)
        v2 = np.ones(u.size)

    cvis = (l1 * v1 + l2 * v2 * phi)/(l1+l2)
    vis2 = np.real(cvis*cvis.conjugate())
    
    return vis2

# =========================================================================
# =========================================================================

def super_gauss(xs, ys, x0, y0, w):
    ''' Returns an 2D super-Gaussian function
    ------------------------------------------
    Parameters:
    - (xs, ys) : array size
    - (x0, y0) : center of the Super-Gaussian
    - w        : width of the Super-Gaussian 
    ------------------------------------------ '''

    x = np.outer(np.arange(xs), np.ones(ys))-x0
    y = np.outer(np.ones(xs), np.arange(ys))-y0
    dist = np.sqrt(x**2 + y**2)

    gg = np.exp(-(dist/w)**4)
    return gg

# =========================================================================
# =========================================================================

def centroid(image, threshold=0, binarize=0):                        
    ''' ------------------------------------------------------
        simple determination of the centroid of a 2D array
    ------------------------------------------------------ '''

    signal = np.where(image > threshold)
    sy, sx = image.shape[0], image.shape[1] # size of "image"
    bkg_cnt = np.median(image)                                       

    temp = np.zeros((sy, sx))
    if (binarize == 1): temp[signal] = 1.0
    else:               temp[signal] = image[signal]

    profx = 1.0 * temp.sum(axis=0)
    profy = 1.0 * temp.sum(axis=1)
    profx -= np.min(profx)                                           
    profy -= np.min(profy)

    x0 = (profx*np.arange(sx)).sum() / profx.sum()
    y0 = (profy*np.arange(sy)).sum() / profy.sum()

    return (x0, y0)

# =========================================================================
# =========================================================================

def find_psf_center(img, verbose=True, nbit=10):                     
    ''' Name of function self explanatory: locate the center of a PSF.

    ------------------------------------------------------------------
    Uses an iterative method with a window of shrinking size to 
    minimize possible biases (non-uniform background, hot pixels, etc)

    Options:
    - nbit: number of iterations (default 10 is good for 512x512 imgs)
    - verbose: in case you are interested in the convergence
    ------------------------------------------------------------------ '''
    temp = img.copy()
    bckg = np.median(temp)   # background level
    temp -= bckg
    mfilt = medfilt(temp, 3) # median filtered, kernel size = 3
    (sy, sx) = mfilt.shape   # size of "image"
    xc, yc = sx/2, sy/2      # first estimate for psf center

    signal = np.zeros_like(img)
    signal[mfilt > 10] = 1.0

    for it in xrange(nbit):
        sz = sx/2/(1.0+(0.1*sx/2*it/(4*nbit)))
        x0 = np.max([int(0.5 + xc - sz), 0])
        y0 = np.max([int(0.5 + yc - sz), 0])
        x1 = np.min([int(0.5 + xc + sz), sx])
        y1 = np.min([int(0.5 + yc + sz), sy])
                                                                     
        mask = np.zeros_like(img)
        mask[y0:y1, x0:x1] = 1.0
        
        #plt.clf()
        #plt.imshow((mfilt**0.2) * mask)
        #plt.draw()

        profx = (mfilt*mask*signal).sum(axis=0)
        profy = (mfilt*mask*signal).sum(axis=1)
        
        xc = (profx*np.arange(sx)).sum() / profx.sum()
        yc = (profy*np.arange(sy)).sum() / profy.sum()
                  
        #pdb.set_trace()
                                                   
        if verbose:
            print("it #%2d center = (%.2f, %.2f)" % (it+1, xc, yc))
            
    return (xc, yc)                                                  

# =========================================================================
# =========================================================================

def recenter(im0, sg_rad=25.0, verbose=True, nbit=10):
    ''' ------------------------------------------------------------
         The ultimate image centering algorithm... eventually...

        im0:    of course, the array to be analyzed
        sg_rad: super-Gaussian mask radius
        bflag:  if passed as an argument, a "bad" boolean is returned
        ------------------------------------------------------------ '''

    szh = im0.shape[1] # horiz
    szv = im0.shape[0] # vertic

    temp = np.max(im0.shape) # max dimension of image

    for sz in [64, 128, 256, 512, 1024, 2048]:
        if sz >= temp: break

    dz = sz/2.           # image half-size

    sgmask = super_gauss(sz, sz, dz, dz, sg_rad)
    x,y = np.meshgrid(np.arange(sz)-dz, np.arange(sz)-dz)
    wedge_x, wedge_y = x*np.pi/dz, y*np.pi/dz
    offset = np.zeros((sz, sz), dtype=complex) # to Fourier-center array

    # insert image in zero-padded array (dim. power of two)
    im = np.zeros((sz, sz))
    orih, oriv = (sz-szh)/2, (sz-szv)/2
    im[oriv:oriv+szv,orih:orih+szh] = im0

    (x0, y0) = find_psf_center(im, verbose, nbit)
    
    im -= np.median(im)

    temp = im * sgmask
    mynorm = temp.sum()

    dx, dy = (x0-dz), (y0-dz)
    im = np.roll(np.roll(im, -int(dx), axis=1), -int(dy), axis=0)

    dx -= np.int(dx)
    dy -= np.int(dy)

    # array for Fourier-translation
    dummy = shift(-dx * wedge_x + dy * wedge_y)
    offset.real, offset.imag = np.cos(dummy), np.sin(dummy)
    dummy = np.abs(shift(ifft(offset * fft(shift(im*sgmask)))))

    #dummy = im
    # image masking, and set integral to right value
    dummy *= sgmask

    return (dummy * mynorm / dummy.sum())

# =========================================================================
# =========================================================================

def get_keck_keywords(hdr):
    '''Extract the relevant keyword information from a fits header.

    This version is adapted to handle NIRC2 data. '''
    data = {
        'tel'    : hdr['TELESCOP'],        # telescope
        'pscale' : 10.0,                   # NIRC2 narrow plate scale (mas)
        'fname'  : hdr['FILENAME'],        # original file name
        'odate'  : hdr['DATE-OBS'],        # UTC date of observation
        'otime'  : hdr['UTC'     ],        # UTC time of observation
        'tint'   : hdr['ITIME'   ],        # integration time (sec)
        'coadds' : hdr['COADDS'  ],        # number of coadds
        'RA'     : hdr['RA'      ],        # right ascension (deg)
        'DEC'    : hdr['DEC'     ],        # declination (deg)
        'filter' : hdr['CENWAVE' ] * 1e-6, # central wavelength (meters)
        # P.A. of the frame (deg) (formula from M. Ireland)
        'orient' : 360+hdr['PARANG']+hdr['ROTPOSN']-hdr['EL']-hdr['INSTANGL']
        }
    print "parang = %.2f, rotposn = %.2f, el=%.2f, instangl=%.2f" % \
        (hdr['PARANG'],hdr['ROTPOSN'],hdr['EL'],hdr['INSTANGL'])
    return data

# =========================================================================
# =========================================================================
def get_nic1_keywords(hdr):
    '''Extract the relevant keyword information from a fits header.

    This version is adapted to handle NICMOS1 data. '''
    data = {
        'tel'    : hdr['TELESCOP'],         # telescope
        'pscale' : 43.1,                    # HST NIC1 plate scale (mas)
        'fname'  : hdr['FILENAME'],         # original file name
        'odate'  : hdr['DATE-OBS'],         # UTC date of observation
        'otime'  : hdr['TIME-OBS'],         # UTC time of observation
        'tint'   : hdr['EXPTIME' ],         # integration time (sec)
        'coadds' : 1,                       # as far as I can tell...
        'RA'     : hdr['RA_TARG' ],         # right ascension (deg)
        'DEC'    : hdr['DEC_TARG'],         # declination (deg)
        'filter' : hdr['PHOTPLAM'] * 1e-10, # central wavelength (meters)
        'orient' : hdr['ORIENTAT'] # P.A. of image y axis (deg e. of n.)
        }
    return data

# =========================================================================
# =========================================================================
def get_idl_keywords(filename):
    '''Extract the relevant keyword information from an idlvar file.
    '''
    data = readsav(filename)
    wavel,bwidth = data['filter']
    data['filter'] = wavel
    data['bwidth'] = bwidth
    data['pscale'] = rad2mas(data.rad_pixel)
    data['fname'] = filename
    data['orient'] = 0
    data['tel'] = 'NIC1'

    return data


# =========================================================================
# =========================================================================
def get_pharo_keywords(hdr):
    '''Extract the relevant keyword information from a fits header.

    This version is adapted to handle PHARO data. '''

    data = {
        'tel'      : hdr['TELESCOP'],         # telescope
        'pscale'   : 25.2,                    # HST NIC1 plate scale (mas)
        'odate'    : hdr['DATE-OBS'],         # UTC date of observation
        'otime'    : hdr['TIME-OBS'],         # UTC time of observation
        'tint'     : hdr['T_INT' ],           # integration time (sec)
        'coadds'   : 1,                       # as far as I can tell...
        'RA'       : hdr['CRVAL1'],           # right ascension (deg)
        'DEC'      : hdr['CRVAL2'],           # declination (deg)
        'filter'   : np.nan, # place-holder   # central wavelength (meters)
        'filtname' : hdr['FILTER'],           # Filter name
        'grism'    : hdr['GRISM'],            # additional filter/nd
        'pupil'    : hdr['LYOT'],             # Lyot-pupil wheel position
        'orient'   : hdr['CR_ANGLE']          # Cassegrain ring angle
        }

    if 'H'       in data['filtname'] : data['filter'] = 1.635e-6
    if 'K'       in data['filtname'] : data['filter'] = 2.196e-6
    if 'CH4_S'   in data['filtname'] : data['filter'] = 1.570e-6
    if 'K_short' in data['filtname'] : data['filter'] = 2.145e-6
    if 'BrG'     in data['filtname'] : data['filter'] = 2.180e-6
    if 'FeII'     in data['grism']   : data['filter'] = 1.648e-6
    
    if np.isnan(data['filter']):
        print("Filter configuration un-recognized. Analysis will fail.")
    return data

# =========================================================================
# =========================================================================
def get_simu_keywords(hdr):
    '''Extract the relevant keyword information from a fits header.

    This is a special version for simulated data. '''
    data = {
        'tel'    : hdr['TELESCOP'],        # telescope
        'pscale' : 11.5,                   # simulation plate scale (mas)
        'fname'  : "simulation",           # original file name
        'odate'  : "Jan 1, 2000",          # UTC date of observation
        'otime'  : "0:00:00.00",           # UTC time of observation
        'tint'   : 1.0,                    # integration time (sec)
        'coadds' : 1,                      # number of coadds
        'RA'     : 0.000,                  # right ascension (deg)
        'DEC'    : 0.000,                  # declination (deg)
        'filter' : 1.6* 1e-6,              # central wavelength (meters)
        'orient' : 0.0                     # P.A. of the frame (deg)
        }
    return data

# =========================================================================
# =========================================================================
def extract_from_array(array, hdr, kpi, save_im=True, re_center=True,
                       wrad=25.0, plotim=False, plotuv=False, wfs=False):
    ''' Extract the Kernel-phase signal from a ndarray + header info.
    
    ----------------------------------------------------------------
    Assumes that the array has been cleaned.
    In order to be able to extract information at the right place,
    a header must be provided as additional argument. This function
    replaces extract_from_fits_frame() when working with a fits
    datacube (multiple frames, one single header).
    
    In addition to the actual Kernel-phase signal, some information
    is extracted from the fits header to help with the interpretation
    of the data to follow.
    
    Parameters are:
    - array: the frame to be examined
    - kpi: the k-phase info structure to decode the data
    - save_im: optional flag to set to False to forget the images
    and save some RAM space
    - wrad: window radius (default = 25 pixels)

    Options:
    -re_center: re-centers the frame before extraction
    - plotim:   plots image
    - plotuv:   plots uv phase map
    - wfs:      wavefront sensing. Instead of kernel-phase, returns phases.

    The function returns a tuple:
    - (kpd_info, kpd_signal)
    - (kpd_info, kpd_signal, im, ac)
    - (kpd_info, kpd_phase)

    ---------------------------------------------------------------- '''

    if 'Keck II' in hdr['TELESCOP']: kpd_info = get_keck_keywords(hdr)
    if 'HST'     in hdr['TELESCOP']: kpd_info = get_nic1_keywords(hdr)
    if 'simu'    in hdr['TELESCOP']: kpd_info = get_simu_keywords(hdr)
    if 'Hale'    in hdr['TELESCOP']: kpd_info = get_pharo_keywords(hdr)
    
    # read and fine-center the frame
    if re_center: im = recenter(array, sg_rad=wrad, verbose=False, nbit=20)
    else:         im = array.copy()

    sz, dz = im.shape[0], im.shape[0]/2  # image is now square

    # meter to pixel conversion factor
    m2pix = mas2rad(kpd_info['pscale']) * sz / kpd_info['filter']

    # rotation of samples according to header info
    th = 90.0 * np.pi/180.
    rmat = np.matrix([[np.cos(th), np.sin(th)], [np.sin(th), -np.cos(th)]])
    uv_rot = np.dot(rmat, kpi.uv.T).T

    uv_samp = kpi.uv * m2pix + dz # uv sample coordinates in pixels
    #uv_samp = uv_rot * m2pix + dz # uv sample coordinates in pixels
    
    # calculate and normalize Fourier Transform
    ac = shift(fft(shift(im)))
    ac /= (np.abs(ac)).max() / kpi.nbh

    xx = np.cast['int'](np.round(uv_samp[:,1]))
    yy = np.cast['int'](np.round(uv_samp[:,0]))
    data_cplx = ac[xx, yy]
    
    vis = np.real(ac*ac.conjugate())
    viscen = vis.shape[0]/2
    vis2 = np.real(data_cplx*data_cplx.conjugate())
    vis2 /= vis[viscen,viscen] #normalise to the origin

    # ---------------------------
    # calculate the Kernel-phases
    # ---------------------------
    kpd_phase = kpi.RED * np.angle(data_cplx) # in radians for WFS
    kpd_signal = np.dot(kpi.KerPhi, kpi.RED*np.angle(data_cplx)) / dtor

    if (save_im): res = (kpd_info, kpd_signal,vis2, im, ac)
    else:         res = (kpd_info, kpd_signal,vis2)
    if (wfs):     res = (kpd_info, kpd_phase)

    uvw = np.max(uv_samp)/2

    if plotim or plotuv:
        plt.clf()
        f0 = plt.subplot(121)
        f0.imshow(im[dz-wrad:dz+wrad,dz-wrad:dz+wrad]**0.5)
        f1 = plt.subplot(122)
        f1.imshow(np.angle(ac))
        f1.plot(uv_samp[:,0], uv_samp[:,1], 'b.')
        f1.axis((dz-uvw, dz+uvw, dz-uvw, dz+uvw))
        plt.draw()
    return res

# =========================================================================
# =========================================================================
def extract_from_fits_frame(fname, kpi, save_im=True, wfs=False):
    ''' Extract the Kernel-phase signal from a single fits frame.
    
    ----------------------------------------------------------------
    Assumes that the fits file has been somewhat massaged, and is
    pretty much ready to go: no pairwise subtraction necessary, 
    etc...
    
    In addition to the actual Kernel-phase signal, some information
    is extracted from the fits header to help with the interpretation
    of the data to follow (think orientation of the telescope!).

    Parameters are:
    - fname: the frame to be examined
    - kpi: the k-phase info structure to decode the data

    Options:
    - save_im: optional flag to set to False to forget the images
    and save some RAM space
    - wfs: wavefront sensing. Instead of kernel-phase, returns phases.

    The function returns a tuple:
    - (kpd_info, kpd_signal)
    - (kpd_info, kpd_signal, im, ac)
    - (kpd_info, kpd_phase)
    ----------------------------------------------------------------  '''

    hdr = pf.getheader(fname)
    if 'Keck II' in hdr['TELESCOP']: kpd_info = get_keck_keywords(hdr)
    if 'HST'     in hdr['TELESCOP']: kpd_info = get_nic1_keywords(hdr)
    if 'simu'    in hdr['TELESCOP']: kpd_info = get_simu_keywords(hdr)
    if 'Hale'    in hdr['TELESCOP']: kpd_info = get_pharo_keywords(hdr)
    
    rev = -1.0
    if 'Hale' in hdr['TELESCOP']: # P3K PA are clockwise
        rev = 1.0
    # read and fine-center the frame
    im = recenter(pf.getdata(fname), sg_rad=40, verbose=False, nbit=40)

    sz = im.shape[0] # image is now square
    dz = sz/2.
    
    # meter to pixel conversion factor
    m2pix = mas2rad(kpd_info['pscale']) * sz / kpd_info['filter']
    uv_samp = kpi.uv * m2pix + dz # uv sample coordinates in pixels

    # calculate and normalize Fourier Transform
    ac = shift(fft(shift(im)))
    ac /= (np.abs(ac)).max() / kpi.nbh

    xx = np.cast['int'](np.round(uv_samp[:,1]))
    yy = np.cast['int'](rev * np.round(uv_samp[:,0]))
    data_cplx = ac[xx, yy]

    vis = np.real(ac*ac.conjugate())
    viscen = vis.shape[0]/2
    vis2 = np.real(data_cplx*data_cplx.conjugate())
    vis2 /= vis[viscen,viscen] #normalise to the origin
    
    # ---------------------------
    # calculate the Kernel-phases
    # ---------------------------
    kpd_phase = kpi.RED * np.angle(data_cplx) # in radians for WFS
    kpd_signal = np.dot(kpi.KerPhi, kpi.RED*np.angle(data_cplx)) / dtor

    if (save_im): res = (kpd_info, kpd_signal, vis2, im, ac)
    else:         res = (kpd_info, kpd_signal, vis2)
    if (wfs):     res = (kpd_info, kpd_phase)

    plotim=True
    wrad = 40.

    uvw = np.max(uv_samp)/2

    if plotim or plotuv:
        plt.clf()
        plt.figure(1, (10,5))
        f0 = plt.subplot(121)
        f0.imshow(im[dz-wrad:dz+wrad,dz-wrad:dz+wrad]**0.5)
        f1 = plt.subplot(122)
        f1.imshow(np.angle(ac))
        f1.plot(uv_samp[:,0], uv_samp[:,1], 'b.')
        f1.axis((dz-uvw, dz+uvw, dz-uvw, dz+uvw))
        plt.draw()
    return res

