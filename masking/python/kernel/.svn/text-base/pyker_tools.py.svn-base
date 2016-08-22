import numpy as np
import pdb

''' ================================================================
    small tools and functions useful for the manipulation of
    Ker-phase data.
    ================================================================ '''

shift = np.fft.fftshift
fft   = np.fft.fft2
ifft  = np.fft.ifft2

''' ---------------------------------------------------------------- 
    convenient little function to convert milliarcsec to radians
    ---------------------------------------------------------------- '''
def mas2rad(x):
    return x*np.pi/(180*3600*1000)

# =========================================================================
# =========================================================================

''' -------------------------------------------------------------------
      Calculate the phases observed by an array on a binary star
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
def phase_binary(u, v, wavel, p):
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

def super_gauss(xs, ys, x0, y0, w):
    ''' -----------------------------------------------------------------
        Returns an array of size (xs,ys), filled with a super-Gaussian
        function centered on (x0,y0) of width (w)
        ----------------------------------------------------------------- '''
    x = np.outer(np.arange(xs), np.ones(ys))-x0
    y = np.outer(np.ones(xs), np.arange(ys))-y0
    dist = np.sqrt(x**2 + y**2)

    gg = np.exp(-(dist/w)**4)
    return gg

# =========================================================================
# =========================================================================

def recenter(im0, th_rad=[80, 40, 10], sg_rad=25.0):
    ''' -----------------------------------------------------------------
        Good image centroid algorithm. Amazing how an apparently simple
        problem turns out to be a real bitch...

        im0:    of course, the array to be analyzed
        th_rad: array of top-hat mask radii to be used
        sg_rad: super-Gaussian mask radius
        ----------------------------------------------------------------- '''

    szh = im0.shape[1] # horiz
    szv = im0.shape[0] # vertic

    temp = np.max(im0.shape) # max dimension of image

    for sz in [64, 128, 256, 512, 1024]:
        if sz > temp: break

    dz = sz/2.           # image half-size
    nm = np.size(th_rad) # number of masks

    sgmask = super_gauss(sz, sz, dz, dz, sg_rad)
    x,y = np.meshgrid(np.arange(sz)-dz, np.arange(sz)-dz)
    dist = np.hypot(y,x)

    wedge_x, wedge_y = x*np.pi/dz, y*np.pi/dz
    offset = np.zeros((sz, sz), dtype=complex) # to Fourier-center array

    # insert image in zero-padded array (dim. power of two)
    im = np.zeros((sz, sz)) 
    orih, oriv = (sz-szh)/2, (sz-szv)/2
    im[oriv:oriv+szv,orih:orih+szh] = im0

    # first round
    # -----------
    profx = im.mean(1) - np.median(im, 1)
    profy = im.mean(0) - np.median(im, 0)
    dx = np.where(profx == profx.max())[0][0] - dz
    dy = np.where(profy == profy.max())[0][0] - dz
    im = np.roll(np.roll(im, -int(dx), axis=0), -int(dy), axis=1)
    print "centroid 0: (dx,dy) = %d %d" % (dx, dy)

    # following rounds: using shrinking mask
    # --------------------------------------
    for i in range(nm):
        thm = np.zeros((sz, sz))               # top-hat mask array
        thm[np.where(dist < th_rad[i])] = 1.0  # define active region
        im1 = im * thm                         # temporary image array
        # ----
        dx = np.sum(np.arange(sz) * im1.sum(1)) / im1.sum() - dz
        dy = np.sum(np.arange(sz) * im1.sum(0)) / im1.sum() - dz
        dx, dy = np.int(np.round(dx)), np.int(np.round(dy))
        im = np.roll(np.roll(im, -dx, axis=0), -dy, axis=1)
        print "centroid %d: (dx,dy) = %d %d" % (i+1, dx, dy)

    # final round: sub-pixel (Fourier) centering
    # ------------------------------------------
    threshold = 0.1

    # create a weighted mask
    wgmask  = np.zeros((sz, sz)) # weight
    wgmask[np.where(im*thm > threshold)] = 1.0
    wgmask *= im

    # determine the centroid offsets
    dx = np.sum(np.arange(sz) * wgmask.sum(1)) / wgmask.sum() - dz
    dy = np.sum(np.arange(sz) * wgmask.sum(0)) / wgmask.sum() - dz
    print "centroid N: (dx,dy) = %+0.3f %+0.3f" % (dx, dy)

    # window function to determine image integral in region of interest
    temp = im * sgmask
    mynorm = temp.sum()

    # array for Fourier-translation
    dummy = shift(dx * wedge_x + dy * wedge_y)
    offset.real, offset.imag = np.cos(dummy), np.sin(dummy)
    dummy = np.abs(shift(ifft(offset * fft(shift(im)))))

    # image masking, and set integral to right value
    dummy *= sgmask
    return (dummy * mynorm / dummy.sum())

# =========================================================================
# =========================================================================

def get_keck_keywords(hdr):
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
    return data

# =========================================================================
# =========================================================================

def get_nic1_keywords(hdr):
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

def load_chi2(fname):
    '''Read in a pickled chi2 grid and associated parameters'''
    import pickle
    
    myf = open(fname,'r')
    data = pickle.load(myf)
    myf.close()
    
    chi2 = data['chi2']
    sepmin = data['sepmin']
    sepmax = data['sepmax']
    anglemin=data['anglemin']
    anglemax = data['anglemax']
    cmin = data['cmin']
    cmax = data['cmax']
    nsep = data['nsep']
    nangle = data['nangle']
    nc = data['nc']
    try:
        null = data['null']
        error = data['error']
    except: null, error = None, None

    return chi2,null,error,sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep,nangle,nc

# =========================================================================
# =========================================================================

def save_chi2(structure,fname,sepmin,sepmax,
              anglemin,anglemax,cmin,cmax,nsep,nangle,nc):
    '''Pickle a chi2 grid and associated parameters'''

    import time, pickle
    
    tic = time.time()
    chi2,null = structure.chi2_3d(sepmin,sepmax,anglemin,anglemax,cmin,cmax,nsep=nsep,
                     nangle = nangle, nc=nc)
    toc = time.time()
    print "Time elapsed =",(toc-tic),"s"

    myf = open(fname,'w')

    data = {
        'chi2'  : chi2,
        'null'  : null,
        'error' : structure.kp_error[0],
        'sepmin': sepmin,
        'sepmax': sepmax,
        'anglemin': anglemin,
        'anglemax': anglemax,
        'cmin'  : cmin,
        'cmax'  : cmax,
        'nsep'  : nsep,
        'nangle': nangle,
        'nc'    : nc}

    pickle.dump(data,myf)
    myf.close()

# =========================================================================
# =========================================================================

def correct_error(structure,chi2):
    '''Correct the errors so chi2/dof comes out to be 1'''

    dof = structure.nkphi - 3

    correction = chi2.min()/dof

    chi2 = chi2/correction

    structure.kp_error = structure.kp_error*np.sqrt(correction)


# =========================================================================
# =========================================================================

def fully_marginalise(chi2,null,error,sepmin,sepmax,anglemin,anglemax,
                      cmin,cmax,nsep,nangle,nc):
    '''Do all the marginalizations to get uncertainties'''

    #first set up the vectors and spacings
    angles = anglemin+(anglemax-anglemin)*np.arange(nangle)/nangle
    seps = sepmin+(sepmax-sepmin)*np.arange(nsep)/nsep
    cons = cmin + (cmax-cmin)*np.arange(nc)/nc
    dxang = (anglemax-anglemin)/nangle
    dxsep = (sepmax-sepmin)/nsep
    dxc = (cmax-cmin)/nc

    #calculate likelihood function
    likelihood = np.exp(-chi2/2.)/np.product(error*np.sqrt(2*np.pi))
    nulllikelihood = np.exp(-null/2.)/np.product(error*np.sqrt(2*np.pi))

    #marginalise once
    margsep = np.trapz(likelihood,x=seps,axis=0)
    margangle = np.trapz(likelihood,x=angles,axis=1)
    margcon= np.trapz(likelihood,x=cons,axis=2)
    
    #marginalise twice
    lsep = np.trapz(margangle,x=cons,axis=1)
    lang = np.trapz(margsep,x=cons,axis=1)
    lc = np.trapz(margangle,x=seps,axis=0)
    
    #find best separation
    bestsep = seps[np.where(lsep==lsep.max())]
    bestang = angles[np.where(lang==lang.max())]
    bestcon = cons[np.where(lc==lc.max())]

    #find normalisation factors
    sepnorm = np.trapz(lsep,x=seps)
    angnorm= np.trapz(lang,x=angles)
    cnorm = np.trapz(lc,x=cons)

    #calculate means

    sepmean = np.trapz(seps*lsep/sepnorm,x=seps)
    angmean = np.trapz(angles*lang/angnorm,x=angles)
    cmean = np.trapz(cons*lc/cnorm,x=cons)
    
    dsep = np.sqrt(np.trapz(seps**2 * lsep/sepnorm,x=seps) - sepmean**2)
    dang = np.sqrt(np.trapz(angles**2 * lang/angnorm,x=angles) - angmean**2)
    dc = np.sqrt(np.trapz(cons**2 * lc/cnorm,x=cons) - cmean**2)

    print 'Separation =',bestsep,'+-',dsep
    print 'Position angle =',bestang,'+-', dang
    print 'Contrast =',bestcon,'+-',dc
    
    return seps,angles,cons,lsep,lang,lc,bestsep,bestang,bestcon,dsep,dang,dc,nulllikelihood
