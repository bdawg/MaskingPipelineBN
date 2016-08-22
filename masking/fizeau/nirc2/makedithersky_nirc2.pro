; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    makedithersky_nirc2.pro                     %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program is intended for use with NIRC2 speckle data.
;     It tries to make up a supersky from a set of dithered data.
;
; Input Variables:
;     dcube    : data cube of frames
; Output:
;    supersky : cleaned, sky data frame
; Options
;    starmask  : size of circular disk to block out speckle splodge
;    darks     : set of skies-supersky
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Written                                         PGT  01Mar04
; Darks added, supersky now median...                       MJI  20Dec06

pro makedithersky_nirc2,dcube,ssky,starmask=starmask,darks=darks,  nocheck = nocheck,  $
 gain = gain, shifts = shifts, handsky=handsky, medsub=medsub
 
;______
; Initialization
;______

if(keyword_set(starmask) eq 0) then starmask=220 ;size of disk to mask out starlight
if(keyword_set(nocheck) eq 0) then nocheck = 0
if(keyword_set(handsky) eq 0) then handsky = 0

info=size(dcube)
if(info[0] eq 2) then begin
   print,'ERROR ### only 1 data frame (no dither!)'
   print,'ERROR ### Cannot make up supersky'
   stop
endif else nsky=info(3)

dimx=info(1) & dimy=info(2)
mask=fltarr(dimx,dimy)
accum=mask
ssky=mask
;darks=replicate(mask,info(3)) ; NO! need size of eventual subarray passed.

; %%%%% Make up circular tapered mask
mask(*)=1
d=shift(dist(dimx,dimy),dimx/2,dimy/2)
w=where(d le starmask/2.)
mask(w)=0
swidth = starmask/10.
mask=smooth(mask,swidth)

; %%%%% Loop over all frames
locs =  fltarr(2, nsky)
mincube=min(dcube, dim=3)
if (nsky lt 3) then medcube = mincube $
else if (nsky lt 6) then $
 medcube = median([ [[mincube]], [[dcube]] ], dim=3) $
else if (nsky lt 8) then $
 medcube = median([ [[mincube]], [[mincube]],[[dcube]] ], dim=3) $
else begin
 mincubes = fltarr(dimx,dimy,nsky/2+1)
 for i=0,nsky/2 do mincubes[*,*,i]=mincube
 medcube = median([ [[mincubes]],[[dcube]] ], dim=3)
endelse
for i=0,nsky-1 do begin
  aframe =  dcube[*, *, i]
  if (keyword_set(medsub)) then aframe -= medcube
  if (keyword_set(gain)) then aframe /= gain
  if (keyword_set(shifts)) then begin
   s =  reverse(sort((shifts[0, i]-shifts[0, *])^2 + (shifts[1, i]-shifts[1, *])^2))
   bias =  dcube[*, *, s[0]]
   if (keyword_set(gain)) then bias /= gain
   aframe -= bias
  endif
  ;Find a standard-deviation in a robust way.
  med =  median(aframe)
  dcube_sdev =  median(abs(aframe-med))*1.5
  ;Smooth the image and find the peak.
  temp = (smooth(aframe,nint(swidth),/edge_truncate))[swidth:dimx-swidth, swidth:dimy-swidth]
  mx=max(temp,mxy)
  mxy = array_indices(temp, mxy)
  xloc=mxy[0]+swidth
  yloc=mxy[1]+swidth
  if (nocheck ne 1) then if (mx lt med + 2*dcube_sdev) then begin
    image_cont,  aframe,  /noc
    print,  'This appears to be a sky-only frame. Type ''S'' to treat as a sky frame...'
;    cursor,  x,  y
;    wait,  0.3
;    if (x gt 0) then begin
    input=''
    read, input
    if (input eq 'S' or input eq 's') then begin
     ssky=ssky+dcube(*,*,i)
     accum += 1.0
     locs[*, i] = [dimx/2, dimy/2]
     continue
    endif
  endif
  if(handsky gt 0) then begin  ; do the whole xy speckle find by eye...
    image_cont,  dcube[*, *, i],  /noc, /asp
    print,  'Manual image finding - click on location of speckle cloud'
    cursor,  xloc,  yloc
    wait,  0.3
  endif
  locs[*, i] =  [xloc,  yloc]
  ssky=ssky+dcube(*,*,i) * shift(mask,xloc-dimx/2,yloc-dimy/2)
  accum=accum+shift(mask,xloc-dimx/2,yloc-dimy/2)
endfor

hole=where(accum eq 0)
if(hole[0] ne -1) then begin
  print,'ERROR ### Could not make up complete sky given dithered frames!'
  image_cont,accum,/nocont,/asp
  stop
endif else ssky=ssky/accum

darks = dcube
;Now for each frame, find an example dark that works.
nused =  intarr(nsky)
for i = 0, nsky-1 do begin
 w =  where(sqrt((locs[0, *] - locs[0, i])^2 + (locs[1, *] - locs[1, i])^2) ge starmask/2)
 if (w[0] eq -1) then stop ;Should NEVER get to here because of error checking above
 s =  sort(nused[w])
 darks[*, *, i] = dcube[*, *, w[s[0]]]
 nused[w[s[0]]] += 1
endfor
;print,  nused

print,'Supersky made up of: ',nsky,' Dithered skies'

end
