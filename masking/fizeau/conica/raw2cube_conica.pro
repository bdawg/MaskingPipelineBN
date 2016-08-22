;######## RELEASE CLIP HERE ######## 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    raw2cube_conica.pro                           %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program is intended for use with CONICA speckle data.
;     It does steps beginning with raw data and outputting a cleaned cube
;
;     At present assumption is that it is passed filenumber of a SINGLE DATA CUBE.
;
; Input Variables:
;     speck.  : vector containing the filename(s) of the source frames.
;     gain    : array containing gain (flat field).
;    datadir  : data directory 
; Output:
;    cube     : cleaned, centered data cube
;  headinfo   : structure of info stripped from header
;  fstats     : frame statistics (flux, sky bgr, xy speckle cloud etc)
;    flux     : array of fluxes for each frame for 10 different
;               aperture sizes.
;  dcube      : cleaned, centered sky cube.
;  quad       : returns the quadrant that the star is found in
; Options:
;    /destripe:     Dstripes the image
;     noskyfit:     Do not fit the sky background, just use supersky.
;  /setsquare :     Set this to trim array size
;  saturation_flag: Returns number of pixels (NOT BAD PIXELS!) which are 
;		      are within 25% of TURNOVER SATURATION! (on AVERAGE per frame)
;		      This usually signals a recent NIRC crash.
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Written (based on nirc2 version)                   PGT  12Jul07

pro raw2cube_conica,speck,ssky,datadir=datadir,dcube=dcube, 	$
        gain=gain,bad_pixels=bad_pixels,skies=skies,tsz=tsz,	$ 
				cube=cube,headinfo=hlog,fstats=fstats,flux=flux,	$
				noskyfit=noskyfit,plotall=plotall,setsquare=setsquare,	$
        nocosmic = nocosmic,  prefix=prefix,extn = extn, 	$
        n_sigma = n_sigma,  dewarp = dewarp,  median_sub = median_sub,  $
        updown=updown, discard_ends=discard_ends, quad=quad,	$
        worry_about_edges=worry_about_edges, cav=cav, pause = pause,    $
	horiz_destripe=horiz_destripe, discard_last=discard_last

if (keyword_set(updown) eq 0) then updown=0
if (keyword_set(prefix) eq 0) then prefix =  'NACO_IMG_SCI'
if (keyword_set(extn) eq 0) then extn =  '_DIT.fits'
if (keyword_set(n_sigma) eq 0) then n_sigma =  5.7 ;An arbitrary number, but works well on test data.
cleanplot,/silent    ; Intialize Display Driver
if (keyword_set(datadir) eq 0) then datadir="./"
if (keyword_set(gain) eq 0) then gain=1.0
if (keyword_set(bad_pixels) eq 0) then bad_pixels=0
if (keyword_set(noskyfit) eq 0) then noskyfit=0
if (keyword_set(setsquare) eq 0) then setsquare=256
if (keyword_set(nocosmic) eq 0) then nocosmic = 0
if (keyword_set(dewarp) eq 0) then dewarp = 0
if (keyword_set(median_sub) eq 0) then median_sub = 0
if (keyword_set(worry_about_edges) eq 0) then worry_about_edges=0
if (keyword_set(cav) eq 0) then cav=1
if (keyword_set(pause) eq 0) then pause=0.1

;set to half chipsize
	if(setsquare eq 1) then setsquare=fix((size(gain))(1)/2)
;______
; Initialization
;______

saturation_flag=long(0)
dimx=(size(gain))(1)
dimy=(size(gain))(2)
skytrim=0.9              ;avoid chip edges when computing sky

; %%%%%  Ensure GAIN map has UNITY Gain for BAD Pixels 
; %%%%%   (important since we later divide by gain.)
gain(bad_pixels)=1.0

; %%%%%  Initialize variables 
proj_x=replicate(1,dimy)
proj_y=replicate(1,dimx)
xpeaks=[0] & ypeaks=[0] & tot_flx=[0.0] & pk_flx=[0.0] & sky_bg=[0.0]
meanxpk=0 & meanypk=0
flux=fltarr(10) ; initialize to accept 10-aperture photometry
han=hanning(dimx,dimy)

; %%%%%  Read in speckle data
filename=datadir+prefix+string(speck,format="(I4.4)")+extn
print,'Reading File: ',filename
incube=float(reform(readfits(filename,head)))
info=size(incube)
nframes=info[3]
dimx=info[1]
dimy=info[2]

; Fix obscure Mar2012 bug with full frame dimy randomly jumping between 1024 and 1026
if(dimy eq 1026) then begin
  dimy = 1024
  incube=incube[*,0:1023,*]
endif

if (keyword_set(discard_ends)) then begin
   nframes -= 1
   incube=incube[*,*,1:nframes-1]
   nframes -=1
endif

if (keyword_set(discard_last)) then begin
   nframes -= 1
   incube=incube[*,*,0:nframes-1]
endif

if (cav gt 0) then begin
    nframes = nframes/cav
    incube = incube[*,*,0:nframes*cav-1]
    incube = rebin(incube, dimx, dimy, nframes)
 endif

hlog=freud(head) 	; read header information

newdimxy=setsquare 
cube  = fltarr(newdimxy, newdimxy,nframes)
dcube = fltarr(newdimxy, newdimxy,nframes)

destripe_stats=fltarr(nframes,4)
shifts =  fltarr(2, nframes)


; Decide if we have a single sky or separate sky for src and cal...
sksz=size(ssky) 
if skies eq -2 then begin                ;; case: separate src/cal skies ...
   if(tsz lt 0) then thissky=ssky[*,*,0]
   if(tsz ge 0) then thissky=ssky[*,*,1]
endif else thissky=ssky                   ; case: single ssky

; %%%%%  Begin loop over data frames 
;t0 = systime(1)
for i=0,nframes-1 do begin

; %%%%%  Subtract off Supersky 
  im=incube[*, *, i]-thissky

; %%%%% Before ANYTHING else, we should remove bad pixels etc.
 im=fix_bad_pixels(im,bad=bad_pixels)

; %%%%% Basic attempt to remove horizontal striping in data (TME May09)
; Median values of each horizontal half-row not containing the
; source are subtracted from the data in the row
if (keyword_set(horiz_destripe)) then begin
  stripes=fltarr(dimy)
  left_tot=total(im[0:dimx/2-1,*])
  right_tot=total(im[dimx/2-1:dimx-1,*])
  if (left_tot lt right_tot) then half_im=im[0:dimx/2-1,*] else half_im=im[dimx/2-1:dimx-1,*]
  for k=0,dimy-1 do begin
    stripes[k]=median(half_im[*,k])
  endfor
  for k=0,dimy-1 do begin
    im[*,k]=im[*,k]-stripes[k]
  endfor
endif

; %%%%% Gain correction:
  im=im/gain

;t1 = systime(1)
;print, 'Time for bad pixels: ', t1-t0
;t0=t1
;%%%%%  Remove cosmic rays or residual bad pixels (MJI) %%%%%
 if (nocosmic eq 0) then begin
   ;7-sigma is lots, but remember that we have 65000 pixels.
   im = sigma_filter_nirc2(im,7,n_sigma=n_sigma,/all,/iterate,/mon)
   ;im = sigma_filter_nirc2(im,7,n_sigma=n_sigma,/all,/iterate)
 endif

;t1 = systime(1)
;print, 'Time for sigma_filter: ', t1-t0
;t0=t1

;Quadrant-by-quadrant offset subtraction. (maybe don't call this for CONICA PGT)
 if (keyword_set(median_sub) ne 0) then begin
  im[0:dimx/2-1, 0:dimy/2-1] -=  median(im[0:dimx/2-1, 0:dimy/2-1]) 
  im[dimx/2:dimx-1, 0:dimy/2-1] -=  median(im[dimx/2:dimx-1, 0:dimy/2-1])
  im[0:dimx/2-1, dimy/2:dimy-1] -=  median(im[0:dimx/2-1, dimy/2:dimy-1])
  im[dimx/2:dimx-1, dimy/2:dimy-1] -=  median(im[dimx/2:dimx-1, dimy/2:dimy-1])
 endif

; % MJI had some fancy shifts and cosmic ray stuff in the Nirc2 version here. Not replicated here.
; Also code for blinking pixels. We can add it if needed to CONCIA later.

; %%%%% Optionally determine D.C. offset in each frame 
; %%%%% Use external frame edge outside 140 pix radius
  if (noskyfit eq 0) then begin             ; Here we do the full sky_subtract
     sky_subtract_nirc,im,140./256.*sqrt(dimx*dimy),stats=stats, $
                      trim=skytrim,plotit=plotall,subtr=0
                      ;trim=skytrim,plotit=plotall,bin=5,subtr=0
  endif else begin
  if (noskyfit eq 1) then begin  ; Here we just calculate the stats (No Subtraction)
     sky_subtract_nirc,im,140./256.*sqrt(dimx*dimy),stats=stats, $
                      trim=skytrim,plotit=plotall,subtr=1
                      ;trim=skytrim,plotit=plotall,bin=5,subtr=1
  endif else stats=[0,0,0,0,0]
  endelse
;t1 = systime(1)
;print, 'Time for sky_subtract_nirc: ', t1-t0
;t0=t1

; %%%%% Determine and Center the speckle cloud. Convolve with a
; Gaussian and find the peak. Don't look near edges.
; (copied from pharomkcube)
bw = 30 ;; border of the detector to be avoided
myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)
temp0 = convol(im, myKer)
if (updown eq 0) then begin
 mx = max(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
 ind = array_indices(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
 xpeak = ind[0]+bw & ypeak = ind[1]+bw
endif else if (updown eq 1) then begin
 mx = max(temp0[bw:dimx-bw, dimy/2+bw:dimy-bw], mxy)
 ind = array_indices(temp0[bw:dimx-bw, dimy/2+bw:dimy-bw], mxy)
 xpeak = ind[0]+bw & ypeak = ind[1]+bw+dimy/2
endif else if (updown eq 2) then begin
 mx = max(temp0[bw:dimx-bw, bw:dimy/2-bw], mxy)
 ind = array_indices(temp0[bw:dimx-bw, bw:dimy/2-bw], mxy)
 xpeak = ind[0]+bw & ypeak = ind[1]+bw
endif else stop
  xpeaks=[xpeaks,xpeak]  ; Recording Peak Position.
  ypeaks=[ypeaks,ypeak]

;t1 = systime(1)
;print, 'Time for peak finding: ', t1-t0
;t0=t1

; %%%%% Check if we are near the chip-edge
  if (worry_about_edges) then if (xpeak gt dimx-newdimxy/2 or xpeak lt newdimxy/2 $
   or ypeak gt dimy-newdimxy/2 or ypeak lt newdimxy/2) then begin
    print,'### WARNING: Center of Speckle Cloud found at :',xpeak,ypeak
    print,'###          SHIFTING - This is too near the edge !!!'
;    im=shift(im,dimx/2-xpeak,dimy/2-ypeak)
;    xpeak=dimx/2 & ypeak=dimy/2
    if (xpeak gt dimx-newdimxy/2) then xpeak=dimx-newdimxy/2
    if (xpeak lt newdimxy/2) then xpeak=newdimxy/2
    if (ypeak gt dimy-newdimxy/2) then ypeak=dimy-newdimxy/2
    if (ypeak lt newdimxy/2) then ypeak=newdimxy/2
  endif

; %%%%% Cut out dark array around periphery
 dkim = (shift(im,-xpeak, -ypeak)) $
  [dimx/2-newdimxy/2:dimx/2+newdimxy/2-1, dimy/2-newdimxy/2:dimy/2+newdimxy/2-1]
; %%%%% Cut data array out 
 im   = (shift(im,dimx/2-xpeak,dimy/2-ypeak)) $
  [dimx/2-newdimxy/2:dimx/2+newdimxy/2-1, dimy/2-newdimxy/2:dimy/2+newdimxy/2-1]

; %%%%% De-Stripe the image (TO BE ADDED?).
; %%%%% Phase 2 De-Stripe of the image removes MULTIPLICATIVE noise
;  stopif (keyword_set(destripe) eq 1) then begin

  if (plotall gt 0) then begin
     if((size(im))[1] gt 256) then tvscl,rebin(im,256,256) else tvscl,im
     legend,['File '+string(speck,format="(I4.4)")+"  Frame "+string(i)]
     ;Give the user a chance to see the image.
     wait,  pause
  endif

; %%%%% Generate Diagnostics 
  tot_flx=[tot_flx,total(im)]
  pk_flx=[pk_flx,max(im[3:newdimxy-4,3:newdimxy-4])] ;;Neglect the edge pixels.
  sky_bg=[sky_bg,stats(1)]

; %%%%% always center image NOT NEEDED AGAIN %%%
;  mx=max(smooth(im,20,/edge),mxy)
;  xpeak=mxy mod dimxx &  ypeak=mxy/dimyy
;  im=shift(im,dimxx/2-xpeak,dimyy/2-ypeak)

; %%%%% Do aperture photometry on centered image
;t1 = systime(1)
;print, 'Time for everything else: ', t1-t0
;t0=t1
  photometry_nirc, im, photometry
  flux=[[flux],[photometry]]

  cube[*,*,i]=im 
  dcube[*,*,i]=dkim
  
;t1 = systime(1)
;print, 'Time for cube-making: ', t1-t0
;t0=t1

endfor                            ;%%%%% FINISH LOOP OVER ALL SPECKLE FRAMES 

; %%%%% Strip useless leading element from Diagnostics:
fstats=fltarr(5,nframes)
fstats(0,*)=xpeaks(1:nframes)
fstats(1,*)=ypeaks(1:nframes)
fstats(2,*)=tot_flx(1:nframes)
fstats(3,*)=pk_flx(1:nframes)
fstats(4,*)=sky_bg(1:nframes)
flux=flux(*,1:nframes)
print,'Mean Location of centroid X,Y: ',mean(fstats[0,*]),', ',mean(fstats[1,*])


;determines the output quadrant
if mean(fstats[0,*]) lt setsquare then $
	if mean(fstats[1,*]) lt setsquare then quad=1 else quad=2 $
	else if mean(fstats[1,*]) lt setsquare then quad=3 else quad=4




; %%%%% Try to find precise parallactic angle
;fix_pa_nirc,headinfo,skyheadinfo,pa,del_pa
;headinfo.pa=pa
;headinfo=create_struct(headinfo,'del_pa',del_pa)

end
