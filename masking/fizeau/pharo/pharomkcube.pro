;; -------------------------------------------------------------------
;; This procedure creates datacubes and their associated cubeinfofile 
;; from a list of pharo file numbers, and other parameters
;;
;; questions: frantz@astro.cornell.edu
;; -------------------------------------------------------------------

pro pharomkcube, filenums, quad, adate, ddir, tsize, savedir, crop = crop, $
                 darkdir=darkdir, flatdir=flatdir, baddir=baddir, $
                 specklecube=specklecube, dither=dither, starmask = starmask, $
                 cal4src=cal4src, firstbad=firstbad, cubedir = cubedir, $
                 suffix=suffix, prefix=prefix, no_fill_names=no_fill_names, $
                 nocosmic = nocosmic, display=display, basic_dither=basic_dither

if NOT keyword_set(cubedir)     then cubedir     = savedir
if NOT keyword_set(flatdir)     then flatdir     = darkdir
if NOT keyword_set(baddir)      then baddir      = darkdir
if NOT keyword_set(cal4src)     then cal4src     = -1
if NOT keyword_set(firstbad)    then firstbad    = 0
if NOT keyword_set(suffix)      then suffix      = '.fits.gz'
if NOT keyword_set(prefix)      then prefix      = 'ph'
if NOT keyword_set(display)     then display     = 0
if NOT keyword_set(dither)      then dither      = 0
if NOT keyword_set(basic_dither)then basic_dither= 0

scsc=uniq(tsize)               ;; This tells how many cubes need to be made
n_cubes = n_elements(scsc)

if(cal4src[0] eq -1) then begin  ;; Default: all cals will be used for each src
   cal4src = intarr(n_cubes, n_cubes) ;; Set this var in flagging later.
   targ  =  where(tsize(scsc) lt 0,  complement = calib)
   if(calib[0] ne -1) then $
     for i =  0, n_elements(targ)-1 do cal4src[calib,targ[i]] =  1  
endif
scsc = [-1,scsc]

frames = filenums
skies = -1
comments='ok'
process_history=adate

;; ---------------------------------------------------------------
;; Total and peak flux rejection criteria. Third parameter rejects
;;frames that are sigma_reject[2] or less of the median peak flux.
;; ---------------------------------------------------------------
sigma_reject =  [4., 4.,  0.4]

;; --------------------------------------------------------
;; Define data structure to hold all the frame statistics:
;; --------------------------------------------------------
ina = intarr(n_cubes)      ;; index of good frames (-1 implies bad frame)
f   = fltarr(n_cubes,2)    ;; this holds frame statistic and err
fl  = fltarr(n_cubes,10,2) ;; average and errors. src#, apert size, 0=mean/1=err

framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}

;; --------------------------------------------------------
;; get the data structure to hold header info and comments
;; and begin to populate it
;; --------------------------------------------------------
olog = make_olog(n_cubes, tsize, frames, skies)

olog.rawdir    = ddir[0]
if keyword_set(cubedir) then olog.cubedir   = cubedir else olog.cubedir = ''
olog.comments  = comments
olog.proc_hist = process_history
olog.cal4src   = cal4src

;; ----------------------------------------------------------------------
;;                start the loop for all data cubes
;; ----------------------------------------------------------------------
for cube_ix = 0, n_cubes-1 do begin

  ;; -------------------------------------------------
  ;; we have to read the header of the first file that
  ;; will go in the data cube
  ;; -------------------------------------------------
  if (n_elements(ddir) eq 1) then scddir = ddir else $
   scddir = ddir[scsc[cube_ix]+1]
  fname = prefix + string(filenums[scsc[cube_ix]+1], $
                          format="(I4.4)") + suffix

  ;; ---------------------
  ;; speckle cube or not ?
  ;; ---------------------
  if (n_elements(specklecube eq 1)) then scspecklecube = specklecube else $
    scspecklecube = specklecube[scsc[cube_ix]+1]

  ;; ---------------------
  ;; how many first bad ?
  ;; ---------------------
  if (n_elements(firstbad eq 1)) then scfirstbad = firstbad else $
    scfirstbad = firstbad[scsc[cube_ix]+1]
  if (quad[0] ge 0) then $
  first = pharoreadfits(scddir+fname, head, /quadrants, $
                        specklecube = scspecklecube, $
                        firstbad = scfirstbad, /silent ) else $
  first = pharoreadfits(scddir+fname, head, specklecube = scspecklecube, $
                        firstbad = scfirstbad, /silent )

  ;; ---------------------
  ;;  find the Flat Field
  ;; ---------------------
  filter = strc(sxpar(head, 'FILTER')) ;; filter to choose the flat
  grism  = strc(sxpar(head, 'GRISM' )) ;; some 'filters' located in grism
  if (n_elements(flatdir) eq 1) then scflatdir = flatdir else $
    scflatdir = flatdir[scsc[cube_ix]+1]
  flatName = filter
  if (filter eq 'CH4_S')       then flatName='H'
  if (filter eq 'Br-gamma')    then flatName='K_short'
  if (filter eq 'CO Bandhead') then flatName='K_short'
  if (filter eq 'Fe II 1.643') then flatName = 'H'
  if (grism  eq 'Fe II 1.643') then flatName = 'H'
  flat = scflatdir + 'flatField.'+flatName+'.fits'

  ;; ---------------------
  ;;  find the badpix map
  ;; ---------------------
  if (n_elements(baddir) eq 1) then scbaddir = baddir else $
    scbaddir = baddir[scsc[cube_ix]+1]
  badpixmap = scbaddir  + 'badpix.fits'

  ;; ---------------------
  ;;  find the badpix map
  ;; ---------------------
  if (n_elements(darkdir) eq 1) then scdarkdir = darkdir else $
    scdarkdir = darkddir[scsc[cube_ix]+1]

  ;; ---------------------
  ;;  dimensions of cube
  ;; ---------------------
  sz = size(first)

  naxis1 = sz[1]
  if not keyword_set(starmask) then starmask = naxis1/2
  naxis2 = sz[2]
  nfiles = scsc[cube_ix+1] - scsc[cube_ix]
  
  if (sz[0] le 3) then depth = 1
  if (sz[0] eq 4) then depth = sz[3]

  stars  = fltarr(naxis1, naxis2, depth * nfiles)
  bckgd  = fltarr(naxis1, naxis2, depth * nfiles)

  ;; ---------------------
  ;; build the cube's name
  ;; ---------------------
   filestring = 'S_'
  if (scspecklecube[0] eq 0) then filestring = 'F_' 
  if (scspecklecube[0] eq 2) then filestring = 'D2_'
  if (scspecklecube[0] eq 4) then filestring = 'D4_'
  filestring += string(filenums[scsc[cube_ix]+1], format="(I4.4)")
  savefile = 'cube' + filestring + '.fits'

  ;; ----------------------------------------
  ;; read information from the fits headers 
  ;; ----------------------------------------
  headinfo = freud(head, olog=olog, ix=cube_ix)  

  if (cube_ix eq 0) then temp1 = olog.mask
  if (olog.mask ne temp1) then begin $
      print,  'ERROR: Multiple masks for the one cubeinfo file not allowed.'
      stop
  endif

  for i = 0, nfiles-1 do begin
      if (n_elements(ddir) eq 1) then scddir = ddir else $
        scddir = ddir[scsc[cube_ix]+1+i]
      fname=scddir + prefix + string(filenums[i + scsc[cube_ix]+1], $
                                     format="(I4.4)") + suffix
      if (quad[0] gt 0) then begin
       a = pharoreadfits(fname, head, /quadrant, specklecube = scspecklecube, $
                        refbias=scdarkdir, badpixf=badpixmap, flatf=flat,  $
                        display=display,  firstbad = scfirstbad, /silent)
       quads = [0,1,2,3]
       theQuad = quad[scsc[cube_ix]+i+1]
       b = where(quads ne theQuad) 

       if (sz[0] eq 3) then begin
           stars[*, *, i*depth:(i+1)*depth-1] = a[*,*,theQuad]
           bg = median(a[*,*,[b[0],b[1],b[2]]], dim=3)
           bckgd[*, *, i*depth:(i+1)*depth-1] = bg
       endif
       if (sz[0] eq 4) then begin
           stars[*, *, i*depth:(i+1)*depth-1] = a[*,*,*,theQuad]
           bg = median(a[*,*,*,[b[0],b[1],b[2]]], dim=4)
           ;stop
           bckgd[*, *, i*depth:(i+1)*depth-1] = bg
           stars[*, *, i*depth:(i+1)*depth-1] -= bg
         endif
       endif else  stars[*, *, i] = $
         pharoreadfits(fname, head, specklecube = scspecklecube, $
                        refbias=scdarkdir, badpixf=badpixmap, flatf=flat,  $
                        display=display,  firstbad = scfirstbad, /silent)
  endfor

  ;; ---------------------
  ;; how about dithering ?
  ;; ---------------------
  if (n_elements(dither eq 1)) then scdither = dither else $
    scdither = dither[scsc[cube_ix]+1]
  if (scdither eq 1) then begin
    makedithersky_nirc2, stars, ssky, starmask=starmask
    for i = 0, nfiles*depth-1 do stars[*,*,i] = stars[*,*,i]-ssky
  endif

  if (n_elements(basic_dither EQ 1 )) THEN scb_dither = basic_dither ELSE $
     scb_dither = basic_dither[scsc[cube_ix]+1]

  if( scb_dither EQ 1 ) THEN BEGIN
     mask_rad = 37
     mask = shift( dist( 2*mask_rad, 2*mask_rad ), mask_rad, mask_rad )
     closed = where( mask LT mask_rad, complement=open )

     mask[ closed ] = !values.f_NAN
     mask[ open   ] = 1.0

     basic_dither, stars, mask, sky, clean_stars, smooth_radius=smooth_radius
     stars = clean_stars
  ENDIF



  ;;----------------------------------------------
  ;; Now time for cosmic ray rejection
  ;; Works well unless there are cosmic rays on the same pixels in
  ;; multiple frames.
  ;;----------------------------------------------
  if not keyword_set(nocosmic) and (nfiles*depth gt 3) then begin
   mncube =  (total(stars,3)-max(stars, dim=3))/(nfiles*depth-1)
   sdevcube =  stars
   for i=0,depth*nfiles-1 do sdevcube[*, *, i] = (stars[*, *, i]-mncube)^2
   sdevcube =  total(sdevcube, 3)-max(sdevcube, dim=3)
   sdevcube = sqrt(sdevcube/(depth*nfiles-2))
   sdevcube = sdevcube > 1.5*smooth(sdevcube,8, /edge)
   for i = 0, nfiles*depth-1 do begin
    w =  where(stars[*, *, i] gt 2*(mncube > 0)+6.*sdevcube)
    if (w[0] ne -1) then begin
      w = array_indices(mncube, w)
      for j = long(0), n_elements(w)/2-1 do $
       stars[w[0, j], w[1, j], i] =  mncube[w[0, j], w[1, j]]
    endif
   endfor
  endif

  ;; -----------------------------------------------------------------
  ;; statistics and photometry of the frames
  ;; fstats: 0 -> xposition of peak
  ;;         1 -> yposition of peak
  ;;         2 -> total frame flux
  ;;         3 -> intensity of peak
  ;;         4 -> sky background
  ;;
  ;; flux: 10-aperture photometry
  ;; -----------------------------------------------------------------
  fstats = fltarr(5, depth * nfiles)
  flux = fltarr(10)

  skymask = fltarr(naxis1, naxis2)
  skymask[*] = 1.0
  xmax = 3*naxis1/4 & xmin = 3*naxis1/4
  ymax = 3*naxis2/4 & ymin = 3*naxis2/4

  for i = xmin, xmax do begin
    for j = ymin, ymax do begin
      skymask[i,j] = 0.0
    endfor
  endfor
  sky = where(skymask gt 0.0)

  ;; -------------------------------------------
  ;; define the gaussian kernel for convol
  ;; -------------------------------------------
  bw = 30 ;; border of the detector to be avoided
  myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)

  for i = 0, nfiles-1 do begin
    for j =  i*depth, i*depth+depth-1 do begin

      ;; --------------------------------
      ;;      center the image
      ;; -------------------------------
      temp0 = convol(stars[*,*,j], myKer)
      mx = max(temp0[bw:naxis1-bw, bw:naxis2-bw], mxy)
      ind = array_indices(temp0[bw:naxis1-bw, bw:naxis2-bw], mxy)
      cx = ind[0]+bw & cy = ind[1]+bw

      stars[*,*,j] = shift(stars[*,*,j],-cx+naxis1/2,-cy+naxis2/2)
      if (j mod 10 eq 0) then print, 'Done: ', j, ' frames.'
      
      ;; -------------------------------
      ;; sky background substraction
      ;; -------------------------------
      temp = stars[*,*,j]
      skylevel = mean(temp[sky])
      print, '[cx,cy] = [',cx,',',cy,'] skylevel = ', skylevel
      stars[*,*,j] = stars[*,*,j] - skylevel

      ;; -----------------------------------
      ;; update the statistics table
      ;; -----------------------------------
      fstats[0,j] = cx
      fstats[1,j] = cy
      fstats[2,j] = total(stars[*,*,j])
      fstats[3,j] = max(stars[*,*,j])
      fstats[4,j] = skylevel
      
      ;; --------------------------------
      ;; photometry 
      ;; -------------------------------
      photometry_nirc, stars[*,*,j], photometry
      flux = [[flux],[photometry]]
    endfor
  endfor
  if keyword_set(crop) then stars =  $
    stars[naxis1/2-crop/2:naxis1/2+crop/2-1,naxis1/2-crop/2:naxis1/2+crop/2-1, *]
  olog.cube_fname[cube_ix,0]=savefile
  olog.cube_fname[cube_ix,1]=filestring
  olog.cube_tsize[cube_ix]=tsize[scsc[cube_ix+1]]

  ;; %%%%% Reduce Aperture Photometry to averages 
  ;; %%%%% and Populate framestats data structure.
  for flx=0,9 do framestats.phot[cube_ix,flx,0] = mean(flux[flx,*])
  for flx=0,9 do framestats.phot[cube_ix,flx,1] = stdev(flux[flx,*])
  framestats.xpk[cube_ix,*]    = [mean(fstats[0,*]), stdev(fstats[0,*])]
  framestats.ypk[cube_ix,*]    = [mean(fstats[1,*]), stdev(fstats[1,*])]
  framestats.totflx[cube_ix,*] = [mean(fstats[2,*]), stdev(fstats[2,*])]
  framestats.pkflx[cube_ix,*]  = [mean(fstats[3,*]), stdev(fstats[3,*])]
  framestats.skybgr[cube_ix,*] = [mean(fstats[4,*]), stdev(fstats[4,*])]

  ;; -------------------------------------------
  ;; Based on the statistics, reject bad frames
  ;; -------------------------------------------
  bad =  where(abs(fstats[2, *]-median(fstats[2, *]))/stdev(fstats[2, *]) gt sigma_reject[0] or $
               abs(fstats[3, *]-median(fstats[3, *]))/stdev(fstats[3, *]) gt sigma_reject[1] or $
               fstats[3, *] lt sigma_reject[2]*median(fstats[3,*]),  complement = good)
  if (bad[0] ne -1) then print,  'Automatically rejecting: ',  n_elements(bad),  ' frames.' else $
    print, 'No bad frames rejected.'
  stars = stars[*, *, good]
  fstats =  fstats[*, good]
  
  olog.cube_sz[cube_ix,*]= (size(stars))[1:3]
  
  ;; ------------------------------
  ;; fill up the frames_ix in olog
  ;; ------------------------------
  olog.frames_ix[0,cube_ix] = frames[scsc[cube_ix]+1]
  olog.frames_ix[1,cube_ix] = frames[scsc[cube_ix+1]]

  ;; -------------------------------------
  ;; write the datacube
  ;; -------------------------------------
  writefits, cubedir+savefile, stars, head
  
endfor

 ;; -------------------------------------------------------------
 ;; Fill in the object names with names from Master_Keck_Catalog
 ;; -------------------------------------------------------------
 if not keyword_set(no_fill_names) then $
   fill_source_name, olog

;; -------------------------------------
;; write the cube info file
;; !!! NB framestats is probably not enough to save here, as we
;; want a record of where each frame's center is. Simply a mean and
;; standard deviation isn't enough for dithered data...
;; -------------------------------------
cinfo = savedir+'cubeinfo'
if (scspecklecube[0] eq 0) then cinfo += '_F_' $
else if (scspecklecube[0] eq 2) then cinfo += '_D2_' $
else if (scspecklecube[0] eq 4) then cinfo += '_D4_' $
else cinfo += '_S_'

save,olog,framestats,file=cinfo + adate + '.idlvar'

end
