; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    raw2cube_nirc2.pro                           %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program is intended for use with NIRC2 speckle data.
;     It does steps beginning with raw data and outputting a cleaned cube
;     Assumes each frame is 2-dimensional (No data cubes)
;
; Input Variables:
;     speck.  : vector containing the filename(s) of the source frames.
;     ssky    : Supersky
;     gain    : array containing gain (flat field).
;    datadir  : data directory 
; Output:
;    cube     : cleaned, centered data cube
;  headinfo   : structure of info stripped from header
;  fstats     : frame statistics (flux, sky bgr, xy speckle cloud etc)
;    flux     : array of fluxes for each frame for 10 different aperture sizes.
; Options:
;    /destripe:     Dstripes the image
;     noskyfit:     Do not fit the sky background, just use supersky.
;  /setsquare :     Set this to trim non-square arrays back to square
;  saturation_flag: Returns number of pixels (NOT BAD PIXELS!) which are 
;		      are within 25% of TURNOVER SATURATION! (on AVERAGE per frame)
;		      This usually signals a recent NIRC crash.
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Written (based on get_full_pow_nxn)                   PGT  28Feb04
;           Extensive changes/cleaning up. 


pro raw2cube_nirc2,speck,ssky,datadir=datadir,  		$
        gain=gain,bad_pixels=bad_pixels,			$
	cube=cube,headinfo=hlog,fstats=fstats,flux=flux, 	$
	noskyfit=noskyfit,setsquare=setsquare,plotall=plotall
	;destripe=destripe,saturation_flag=saturation_flag

cleanplot,/silent    ; Intialize Display Driver
if (keyword_set(datadir) eq 0) then datadir="./"
if (keyword_set(gain) eq 0) then gain=1.0
if (keyword_set(bad_pixels) eq 0) then bad_pixels=0
if (keyword_set(noskyfit) eq 0) then noskyfit=0

;______
; Initialization
;______

saturation_flag=long(0)
dimx=(size(gain))(1)
dimy=(size(gain))(2)
nframes=n_elements(speck)
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
destripe_stats=fltarr(nframes,4)

; Define data structure to hold Header info 
s=replicate('',nframes) &  i=intarr(nframes) &  f=fltarr(nframes) &  d=dblarr(nframes)
hlog={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,$ ;from headers
        optic_cfg:s,lyot:s,grism:s,source_name:s,utc:s,                 $ 
        date:s,jd:d,elevation:f,airmass:f,pa:f,del_pa:d}

; %%%%%  Begin loop over data frames 
for i=0,nframes-1 do begin

; %%%%%  Read in speckle data
  filename=datadir+"n"+string(speck[i],format="(I4.4)")+".fits"
  print,'Reading File: ',filename
  frame=float(reform(readfits(filename,head)))
  info=size(frame)
  headinfo=freud(head) 	; read header information

; Now Check for 85% Saturated Pixels.
;  number=where( (frame ge saturateval) and (bad_mask ne 1), ct)
;  saturation_flag=saturation_flag+ct
;  if (ct gt 0) then begin
;     print,'Number of Pixels at 85% of total Saturation: ',ct
;  endif

; %%%%%  Subtract off Supersky 
  im=(frame-ssky)

; %%%%% Gain correction:
  im=im/gain

; %%%%%  Remove Bad Pixels (Preliminary Pass)
  im=fix_bad_pixels(im,bad=bad_pixels)

  
; %%%%% Populate hlog data structure.
  hlog.instrument[i] =headinfo.instrument
  hlog.nax1[i]       =headinfo.nax1
  hlog.nax2[i]       =headinfo.nax2
  hlog.t_int[i]      =headinfo.t_int
  hlog.coadd[i]      =headinfo.coadd
  hlog.filter[i]     =headinfo.filter
  hlog.slit[i]       =headinfo.slit
  hlog.optic_cfg[i]  =headinfo.optic_cfg
  hlog.lyot[i]       =headinfo.lyot
  hlog.grism[i]      =headinfo.grism
  hlog.source_name[i]=headinfo.source_name
  hlog.utc[i]        =headinfo.utc
  hlog.date[i]       =headinfo.date
  hlog.jd[i]         =headinfo.jd
  hlog.elevation[i]  =headinfo.elevation
  hlog.airmass[i]    =headinfo.airmass
  hlog.pa[i]         =headinfo.pa


; %%%%% Optionally determine D.C. offset in each frame 
; %%%%% Use external frame edge outside 140 pix radius
  if (noskyfit eq 0) then begin             ; Here we do the full sky_subtract
     sky_subtract_nirc,im,140./256.*sqrt(dimx*dimy),stats=stats, $
                      trim=skytrim,plotit=plotall,bin=5,subtr=0
  endif else begin
  if (noskyfit eq 1) then begin  ; Here we just calculate the stats (No Subtraction)
     sky_subtract_nirc,im,140./256.*sqrt(dimx*dimy),stats=stats, $
                      trim=skytrim,plotit=plotall,bin=5,subtr=1
  endif else stats=[0,0,0,0,0]
  endelse
  if (plotall gt 0) then begin
     tvscl,im
     legend,['File n'+string(speck[i],format="(I4.4)")+".fits"]
  endif

; %%%%% Determine and Center the speckle cloud (use Center of Gravity)
  profile=total(im,2)
  xpeak=total(profile*findgen(n_elements(profile))) / total(profile)
  profile=total(im,1)
  ypeak=total(profile*findgen(n_elements(profile))) / total(profile)
  if (xpeak gt 200./256*dimx or xpeak lt 50./256*dimx $
    or ypeak gt 200./256*dimy or ypeak lt 50./256*dimy) then begin
    print,'### WARNING: Center of Speckle Cloud found at :',xpeak,ypeak
  endif
  xpeaks=[xpeaks,xpeak]  ; Recording Peak Position.
  ypeaks=[ypeaks,ypeak]


; %%%%% Optionally trim array square
  if (keyword_set(setsquare) eq 1) then begin
     if (setsquare gt 1) then newdimxy=setsquare $
     else newdimxy=min([dimx,dimy])
     im=grabnxn(im,newdimxy,x=dimx/2,y=dimy/2)
     dimx=(size(im))[1] & dimy=(size(im))[2]
     if(dimx ne newdimxy or dimy ne newdimxy) then $
       print,'WARNING### did not square array properly'
  endif

; %%%%% De-Stripe the image (TO BE ADDED?).
; %%%%% Phase 2 De-Stripe of the image removes MULTIPLICATIVE noise
;  if (keyword_set(destripe) eq 1) then begin


; %%%%% Generate Diagnostics 
  tot_flx=[tot_flx,total(im)]
  pk_flx=[pk_flx,max(im)]
  sky_bg=[sky_bg,stats(1)]


; %%%%% center image
  im=shift(im,dimx/2-xpeak,dimy/2-ypeak)

; %%%%% Do aperture photometry on centered image
  photometry_nirc, im, photometry
  flux=[[flux],[photometry]]

; %%%%% Overwrite the data cube with the cleaned frames
  if(i eq 0) then cube=im else cube=[[[cube]],[[im]]]

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
