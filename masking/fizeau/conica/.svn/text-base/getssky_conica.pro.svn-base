;######## RELEASE CLIP HERE ######## 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    getssky_conica.pro                           %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
;
; Input Variables:
;     speck   : vector containing the filename(s) of the source frames.
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
;   new_bad_pixels: (optional output) An array showing the bad pixels
;                   that were continually flagged as cosmic rays, so
;                   that the flat field bad pixel map can be updated
;   discard_ends:   Throw away the first and last frame from every
;                   data cube. The last frame is a sum of the rest,
;                   and the first suffers from reset issues (NACO manual)
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Written (based on nirc2 version)                   PGT  12Jul07

pro getssky_conica,speck,ssky,datadir=datadir,dcube=dcube, 	bad_pixels=bad_pixels,$
        prefix=prefix,extn = extn,gain=gain,starmask=starmask,handsky=handsky, $
        sky_lr=sky_lr, discard_ends=discard_ends, polzflag=polzflag, medsub=medsub,new_bad_pixels=new_bad_pixels

if (keyword_set(prefix) eq 0) then fprefix = replicate('NACO_IMG_SCI', n_elements(speck))
if (keyword_set(extn) eq 0) then extn =  '_DIT.fits'
if (keyword_set(datadir) eq 0) then datadir="./"

; if (n_elements(prefix) lt n_elements(speck)) then fprefix = replicate(prefix, n_elements(speck)) else fprefix=prefix

;______
; Initialization
;______

dimx=(size(gain))(1)
dimy=(size(gain))(2)
nframes=n_elements(speck)

ssky=fltarr(dimx,dimy)

n_files=n_elements(speck)
new_bad_pixels=fltarr(dimx,dimy)
for cube_ix=0,n_files-1 do begin 

; %%%%%  Read in speckle data
  filename=datadir+prefix[cube_ix]+string(speck[cube_ix],format="(I4.4)")+extn
  incube=float(reform(readfits(filename,head)))
; ; Fix obscure Mar2012 bug with full frame dimy randomly jumping between 1024 and 1026
  if( (size(incube))[2] eq 1026) then incube=incube[*,0:1023,*]

  nframes=(size(incube))[3]
  if discard_ends eq 1 then begin
     nframes -= 1
     incube = incube[*,*,1:nframes-1]
     nframes -=1
  endif else if discard_ends eq 2 then begin
     nframes-=1
     incube=incube[*,*,0:nframes-1]
  endif
  incube=total(incube,3)/nframes
  incube=fix_bad_pixels(incube,bad=bad_pixels)
  ;PGT: first get rid of cosmic/badpx_spikes before smoothing
  ;MJI: moved to here from make_dither_sky_nirc2
  incube=sigma_filter_nirc2(incube,7,n_sigma=n_sigma,/all,/iterate,bad_pixels=wbad)
  image_cont,incube,/n,/a

  ;;add the bad pixels from sigma_filter_nirc2 together, so that we
  ;;can pick up any bad pixels missed in the flats. We may want to use
  ;;this to update the flats
  if wbad[0] ne -1 then new_bad_pixels[wbad]+=1

  if(cube_ix eq 0) then bigcube=incube else bigcube=[[[bigcube]],[[incube]]]
endfor

;;normalise the new bad pixel map so it is shown as the % of the
;;frames that each pixel was marked as bad
new_bad_pixels/=n_files

;print,'Click' ;Testing
;cursor,xxxx,yyyy ;Testing

; %%%%% Make up supersky
if (keyword_set(sky_lr) ) then begin
    nx = (size(incube))[1]
    ssky=incube
    ssky[*]=0
    nfr=ssky
    for i=0,n_files-1 do begin
      ;profile = total(bigcube[*,*,i], 2)
      the_median=median(bigcube[*,*,i])
      profile = total(bigcube[*,*,i]-gain*the_median, 2)
      m=max(profile, ix)
      if (ix gt nx/2) then begin
          nfr[0:nx/2-1,*] += 1
          ssky[0:nx/2-1,*] += bigcube[0:nx/2-1,*,i]
          ;image_cont,bigcube[0:nx/2-1,*,i],/nocont,/asp ;Testing
      endif else begin
          nfr[nx/2:*,*] += 1
          ssky[nx/2:*,*] += bigcube[nx/2:*,*,i]
          ;image_cont,bigcube[nx/2:*,*,i],/nocont,/asp ;Testing
          ;print,speck[i] ;Testing
      endelse
      ;wait,1.0 ;Testing
    endfor
  if (min(nfr) eq 0) then stop
  ssky /= nfr
endif else begin
 
;stop
makedithersky_conica,bigcube,ssky,starmask=starmask,polzflag=polzflag,darks=darks,gain=gain,medsub=medsub
endelse

;Write out the ssky to view later... !no! do this later - could be more skies to make.
;save,ssky,file='supersky_saved.idlvar'

end
