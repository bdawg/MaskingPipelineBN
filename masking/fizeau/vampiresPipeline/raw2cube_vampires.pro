;######## RELEASE CLIP HERE ######## 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    raw2cube_vampires.pro                           %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Based on raw2cube_conica.pro. Creates 'cleaned' cubes from VAMPIRES
; data.
; BN Aug2013
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
; Options:
;    /destripe:     Dstripes the image
;     noskyfit:     Do not fit the sky background, just use supersky.
;  /setsquare :     Set this to trim array size
;  saturation_flag: Returns number of pixels (NOT BAD PIXELS!) which are 
;		      are within 25% of TURNOVER SATURATION! (on AVERAGE per frame)
;		      This usually signals a recent NIRC crash.
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro raw2cube_vampires,speck,ssky,datadir=datadir,dcube=dcube, 	$
        gain=gain,bad_pixels=bad_pixels,skies=skies,tsz=tsz,	$
	cube=cube,headinfo=hlog,fstats=fstats,flux=flux, 	$
	noskyfit=noskyfit,plotall=plotall,setsquare=setsquare, 	$
        nocosmic = nocosmic,  prefix=prefix,extn = extn, 	$
        n_sigma = n_sigma,  dewarp = dewarp,  median_sub = median_sub, $
        updown=updown, discard_last=discard_last, $
        worry_about_edges=worry_about_edges, cav=cav,horiz_destripe=horiz_destripe, $
        lcvrstate=lcvrstate, skyrad=skyrad, speckpos=speckpos, keepzeronum=keepzeronum
     
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

;if (keyword_set(discard_last) eq 0) then discard_last=0
if (keyword_set(horiz_destripe) eq 0) then horiz_destripe=0

if (setsquare eq 1) then setsquare=256
;______
; Initialization
;______
common metadatatype

saturation_flag=long(0)
dimx=(size(gain))(1)
dimy=(size(gain))(2)
;skytrim=0.9              ;avoid chip edges when computing sky
skytrim=0

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
if keyword_set(keepzeronum) eq 0 then begin
if speck eq 0 then begin
    fnumstring = ''
endif else begin
    fnumstring = strn(speck)
endelse
endif else begin
    fnumstring = strn(speck)
 endelse


;filename=datadir+prefix+string(speck,format="(I4.4)")+extn
filename=datadir+prefix+fnumstring+extn
print,'Reading File: ',filename
incube=float(reform(readfits(filename,head)))
; For Post-2013 VAMPIRES data - read extension 1 to get proper header data
if header_method eq 1 then dummy=readfits(filename,headExtn1,exten_no=1)

; %%%%%  As per LCVR state, only include appropriate frames
; This assumes the FT (July2013) method, where interleaved frames
;     are discarded.
tnf=(size(incube))(3)
if lcvrstate eq 1 then begin
    p_inds=indgen(tnf/4)*4+1 ;Cube indices for LCVR state 1
endif else begin
    p_inds=indgen(tnf/4)*4+3 ;Cube indices for LCVR state 2
endelse
incube=incube[*,*,p_inds]
; %%%%%

info=size(incube)
nframes=info[3]
dimx=info[1]
dimy=info[2]
if (keyword_set(discard_last)) then begin
   nframes -= 1
   incube=incube[*,*,0:nframes-1]
endif

if (cav gt 0) then begin
    nframes = nframes/cav
    incube = incube[*,*,0:nframes*cav-1]
    incube = rebin(incube, dimx, dimy, nframes)
 endif



;hlog=freud(head) 	; read header information
;;; Different ways to get header info for VAMPIRES...

case header_method of 
   0: begin
      s='' &  i=0 &  f=float(0.0) &  d=double(0.0)
      headinfo={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,optic_cfg:s, $
          lyot:s,grism:s,source_name:s,utc:s,date:s,jd:d,elevation:f,airmass:f,pa:f, $
          ra:d,dec:d,equinox:f, mask:s,  raoff:f, decoff:f,  del_elev:f,  del_pa:f, $
          emgain:f, hwp:f, timingpattern:s, imgRotAng:f, imgRotPad:f, imgRotPap:f ,$
          ADCstagepos:f, ADCp1angle:f, ADCp2angle:f, azimuth:f, localtime:s} ;Last 2 lines are added for VAMPIRES
      hlog=headinfo
      end
   
   1: begin
      print,'Using header method 1...'
      hlog = jung(headExtn1)

      end

   2: begin
      ; Use Lacan to read signifiers
      ; and return them in a structure...
      ; Assumes VAMPIRES log is in datadir/LogfileVampires.txt
      fitsfilename=prefix+fnumstring+extn
      hlog=lacan(datadir,fitsfilename)
      end

endcase

hlog.nax1=dimx
hlog.nax2=dimy
hlog.nax3=info[3]
hlog.instrument[*] = 'VAMPIRES'

; #### TO DO: Populate these things:
; mask = '9hole' ;####
; filter = '750-40' ;####
; pa = 0.0 ;####
; hlog.mask[*] = mask
; hlog.filter[*]=filter
; hlog.pa[*]= pa




newdimxy=setsquare 
cube  = fltarr(newdimxy, newdimxy,nframes)
dcube = fltarr(newdimxy, newdimxy,nframes)

destripe_stats=fltarr(nframes,4)
shifts =  fltarr(2, nframes)

; Decide if we have a single sky or separate sky for src and cal...
sksz=size(ssky) 
if(sksz[0] eq 3) then begin               ; case: separate src/cal skies ...
   if(tsz lt 0) then thissky=ssky[*,*,0]
   if(tsz ge 0) then thissky=ssky[*,*,1]
endif else thissky=ssky                   ; case: single ssky

; %%%%%  Begin loop over data frames 
;t0 = systime(1)
for i=0,nframes-1 do begin

; %%%%%  Subtract off Supersky 
  im=incube[*, *, i]-thissky

; %%%%% Before ANYTHING else, we should remove bad pixels etc.
; im=fix_bad_pixels(im,bad=bad_pixels) ; Saying no bad pix with emccd

; %%%%% Basic attempt to remove horizontal striping in data (TME May09)
; Median values of each horizontal half-row not containing the
; source are subtracted from the data in the row
; if (keyword_set(horiz_destripe)) then begin
;   stripes=fltarr(dimy)
;   left_tot=total(im[0:dimx/2-1,*])
;   right_tot=total(im[dimx/2-1:dimx-1,*])
;   if (left_tot lt right_tot) then half_im=im[0:dimx/2-1,*] else half_im=im[dimx/2-1:dimx-1,*]
;   for k=0,dimy-1 do begin
;     stripes[k]=median(half_im[*,k])
;   endfor
;   for k=0,dimy-1 do begin
;     im[*,k]=im[*,k]-stripes[k]
;   endfor
; endif

; %%%%% Gain correction:
  im=im/gain

;t1 = systime(1)
;print, 'Time for bad pixels: ', t1-t0
;t0=t1

;%%%%%  Remove cosmic rays or residual bad pixels (MJI) %%%%%
 if (nocosmic eq 0) then begin
   ;7-sigma is lots, but remember that we have 65000 pixels.
   im = sigma_filter_nirc2(im,7,n_sigma=n_sigma,/all,/iterate,/mon)
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
; Old conica/nirc used external frame edge outside 140 pix radius... 
; Now radius is provided, and two masks are made.
  if (noskyfit eq 0) then begin             ; Here we do the full sky_subtract
     sky_subtract_vampires,im,skyrad,stats=stats, $
                      trim=skytrim,plotit=plotall,subtr=0
                      ;trim=skytrim,plotit=plotall,bin=5,subtr=0
  endif
  if (noskyfit eq 1) then begin  ; Here we just calculate the stats (No Subtraction)
     sky_subtract_vampires,im,skyrad,stats=stats, $
                      trim=skytrim,plotit=plotall,subtr=1
                      ;trim=skytrim,plotit=plotall,bin=5,subtr=1
  endif
  if (noskyfit eq 200) then begin
      im=im-200.
      print,'Subtracting bias of 200 ADU, no sky subtraction'
      stats=[0,0,0,0,0]
  endif 
  if (noskyfit eq 2) then begin
      stats=[0,0,0,0,0]
  endif
  if (noskyfit eq 3) then begin
     sky_subtract_vampires,im,skyrad,stats=stats, $
         trim=skytrim,plotit=plotall,subtr=0,submean=1
  endif


;t1 = systime(1)
;print, 'Time for sky_subtract_nirc: ', t1-t0
;t0=t1

; %%%%% Determine and Center the speckle cloud. Convolve with a
; Gaussian and find the peak. Don't look near edges.
; (copied from pharomkcube)
; Note to self: Vis speckles are huge and centre is poorly defined. 
; Does it calibrate better if we use the (xpk,ypk) posn for the whole
; cube or even the whole set? 
if speckpos[0] eq -1 then begin
    bw = 30 ;; border of the detector to be avoided
    ;myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)
    ;myKer=psf_gaussian(fwhm=50,npixel=60) ;Big and slow but does the trick...
    ;temp0 = convol(im, myKer)
    temp0=smooth(im,50) ;Quick and dirty
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

endif else begin
    if updown eq 0 then begin
        print,"Updown is 0 - this doesn't seem right..."
    endif
    if updown eq 1 then begin
        xpeak = speckpos[0,0]
        ypeak = speckpos[1,0]
    endif
    if updown eq 2 then begin
        xpeak = speckpos[0,1]
        ypeak = speckpos[1,1]
    endif
endelse

xpeaks=[xpeaks,xpeak]       ; Recording Peak Position.
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

; Hack to show speckle well in nx
;wset,1 ;;;
  if (plotall gt 0) then begin
     if((size(im))[1] gt 256) then tvscl,rebin(im,256,256) else tvscl,im
     legend,['File '+string(speck,format="(I4.4)")+"  Frame "+string(i)]
     ;Give the user a chance to see the image.
     wait,  0.001 ; No, I'm in a hurry!
  endif
;wset,0 ;;;

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

; %%%%% Try to find precise parallactic angle
;fix_pa_nirc,headinfo,skyheadinfo,pa,del_pa
;headinfo.pa=pa
;headinfo=create_struct(headinfo,'del_pa',del_pa)

end
