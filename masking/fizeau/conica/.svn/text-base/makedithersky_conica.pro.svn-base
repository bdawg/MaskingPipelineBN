;######## RELEASE CLIP HERE ######## 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;      makedithersky_conica
;
; PURPOSE:
;      Construct a supersky from a set of dithered polarised Conica
;      data.
;
; EXPLANATION:
;      Combines the Numerical Recipes functions SPL_INIT and SPL_INTERP
;
; CALLING SEQUENCE:
;      makedithersky_conica,bigcube,ssky,starmask=starmask,darks=darks
;
; INPUTS:
;      bigcube - data cube constructed from the average counts of the raw data cubes.  
;
; INPUT-OUTPUT KEYWORD:
;      starmask - diameter of the circular mask used to blank out the source.
;      darks - 
;
; OUTPUT:
;       ssky - supersky made up of the average sky counts for each pixel. 
;
; HISTORY:
;      version 1.0  TME 28May09  Adapted from makedithersky_nirc2.pro
;      version 2.0  MJI 13Apr10  Added the NIRC2 features. THIS CODE
;      HAS DIVERGED AND IS STUPID!!! TODO: RE-CONVERGE CONICA AND
;      NIRC2 CODES!!!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



pro makedithersky_conica,dcube,ssky,starmask=starmask,polzflag=polzflag,darks=darks,gain=gain, medsub=medsub
 
;______
; Initialization
;______

if(keyword_set(starmask) eq 0) then starmask=230 ;size of disk to mask out starlight
if(keyword_set(polzflag) eq 0) then polzflag=0   ; 1=search-and-blank 2 stars for each dither 

info=size(dcube)
if(info[0] eq 2) then begin
   print,'ERROR ### only 1 data frame (no dither!)'
   print,'ERROR ### Cannot make up supersky'
   stop
endif else nsky=info(3)

dimx=info(1) & dimy=info(2)
starblock=fltarr(dimx*2,dimy*2)
ssky=fltarr(dimx,dimy)
accum=ssky
;darks=replicate(mask,info(3)) ; NO! need size of eventual subarray passed.

;; Make up circular tapered mask
starblock(*)=1
d=shift(dist(dimx*2,dimy*2),dimx,dimy)
w=where(d le starmask/2.)
starblock(w)=0
swidth = starmask/10.
starblock=(smooth(starblock,swidth))

;; Loop over all frames
locs1 =  fltarr(2, nsky)
locs2 = fltarr(2, nsky)

mincube=min(dcube, dim=3)
if (nsky lt 3) then medcube = mincube $
else if (nsky lt 6) then $
 medcube = median([ [[mincube]], [[dcube]] ], dim=3) $
else if (nsky lt 8) then $
 medcube = median([ [[mincube]], [[mincube]],[[dcube]] ], dim=3) $
else begin
 mincubes = fltarr(dimx,dimy,nsky/2-1)
 for i=0,nsky/2-2 do mincubes[*,*,i]=mincube
 medcube = median([ [[mincubes]],[[dcube]] ], dim=3)
endelse

for i=0,nsky-1 do begin
  ;if (keyword_set(medsub)) then aframe -= medcube
  aframe =  dcube[*, *, i] 
  if (keyword_set(medsub)) then aframe -= medcube
  if(keyword_set(gain)) then aframe=aframe-gain*median(aframe)
  temp = (smooth(aframe,nint(swidth),/edge))[swidth:dimx-swidth, swidth:dimy-swidth]
  ;; Find the star (or first polarisation) on the chip
  mx1=max(temp,mxy1)
  mxy1 = array_indices(temp, mxy1)
  xloc1=mxy1[0]+swidth
  yloc1=mxy1[1]+swidth       
  locs1[*, i] =  [xloc1,  yloc1]
  if(polzflag eq 0) then begin    ; case for single star (non-wollaston)
      thisblock=(shift(starblock,xloc1-dimx/2,yloc1-dimy/2))[dimx/2:3*dimx/2-1,dimy/2:3*dimy/2-1]
      ssky=ssky+dcube(*,*,i) * thisblock
      accum=accum+thisblock
  endif else begin                ; case for polarized data
  ;; Cut the first polarisation out of the chip; fixed 16042010
      temp = aframe*(shift(starblock,xloc1-dimx/2,yloc1-dimy/2))[dimx/2:3*dimx/2-1,dimy/2:3*dimy/2-1]
      temp = ((smooth(temp,nint(swidth),/edge))[swidth:dimx-swidth, swidth:dimy-swidth])      
  ;; With the first polarisation removed, find the second polarisation on the chip
      mx2=max(temp,mxy2)
      mxy2 = array_indices(temp, mxy2)
      xloc2=mxy2[0]+swidth
      yloc2=mxy2[1]+swidth       
      locs2[*, i] =  [xloc2,  yloc2]
      thisblock=(shift(starblock,xloc1-dimx/2,yloc1-dimy/2))[dimx/2:3*dimx/2-1,dimy/2:3*dimy/2-1] * $
                (shift(starblock,xloc2-dimx/2,yloc2-dimy/2))[dimx/2:3*dimx/2-1,dimy/2:3*dimy/2-1]
      ssky=ssky+dcube(*,*,i) * thisblock
      accum=accum+thisblock
  endelse
endfor

sky_gap=where(accum eq 0)
if(sky_gap[0] ne -1) then begin
  print,'ERROR ### Could not make up complete sky given dithered frames!'
  image_cont,accum,/nocont,/asp
  stop
endif else ssky=ssky/accum

;; darks = dcube
;; ;Now for each frame, find an example dark that works.
;; nused =  intarr(nsky)
;; for i = 0, nsky-1 do begin
;;  w =  where(sqrt((locs[0, *] - locs[0, i])^2 + (locs[1, *] - locs[1, i])^2) ge starmask/2)
;;  if (w[0] eq -1) then stop ;Should NEVER get to here because of error checking above
;;  s =  sort(nused[w])
;;  darks[*, *, i] = dcube[*, *, w[s[0]]]
;;  nused[w[s[0]]] += 1
;; endfor
;; ;print,  nused

print,'Supersky made up of: ',nsky,' Dithered skies'

end
