;+
; NAME:
;      centerconstell.pro
; PURPOSE:
; 	TO find (and correct) global tip-tilt offset for 
;       constellation pattern of LWS seg tilt data.
;
; USAGE:
;   centerconstell
;
; INPUTS:
;	cube    - data cube
;
; RETURNS
;       shifted cube
;
; OPTIONAL INPUTS:
;	pattern	- one of the known Keck ACS segment patterns
;
; OUTPUTS:
;  cents	- Vector of centers
;
; Written					28Aug04 PGT

function centerconstell,cube,pattern=pattern,cents=cents

ccube=cube
pix_scale=82.7              ; compromise plate scale milli-arcseconds/pixel
if (keyword_set(pattern) eq 0) then pattern=0

if(pattern eq 0) then pxy=[0,0]

if(pattern eq 1) then pxy=[[1.97070, -3.41334],[-3.94139, 0.0],[1.97070, 3.41334],[0.0,0.0]]
if(pattern eq 2) then pxy=[0,0]
if(pattern eq 3) then pxy=[[-5.0, -5.0],[0.0, -5.0],[0.0, 0.0],[-5.0,0.0]] + 2.5
if(pattern eq 4) then pxy=[[1.97070, -3.41334],[-3.94139, 0.0],[1.97070, 3.41334],[0.0,0.0]]
if(pattern eq 5) then pxy=[0,0]
if(pattern eq 6) then pxy=[[-5.0, -5.0],[0.0, -5.0],[0.0, 0.0],[-5.0,0.0]] + 2.5

; convert to pixels on chip
sz=size(cube)
pxy=nint(pxy*1000./pix_scale)
pxy=pxy+sz[1]/2

sp=size(pxy)
if(sp[0] eq 1) then n_spots=1 else n_spots=sp[2]
xcorf=fltarr(sz[1],sz[2])


; %%%%% Make up circular gaussian spot
gspot=xcorf
d=dist(sz[1],sz[2])
spotsize=3               ;fwhm of gaussian to make up
g=gaussian(findgen(sz[1]),[1.0,0,spotsize])
gspot = g(round(d))
for p=0,n_spots-1 do xcorf=xcorf+shift(gspot,pxy[0,p],pxy[1,p])
ftxc=conj(fft(xcorf,1))

; %%%%% Find shifts by correlation
cents=intarr(2,sz(3))
for ff=0,sz(3)-1 do begin 
  ftc=fft(cube(*,*,ff),1)
  xc=real_part(fft(ftc*ftxc))
  m=max(shift(xc,sz[1]/2,sz[2]/2),midx)
  x_cent=(midx mod sz[1] ) - sz[1]/2
  y_cent=(midx/sz[2] ) - sz[2]/2
  cents[*,ff]=[x_cent,y_cent]
  ccube(*,*,ff)=shift(cube(*,*,ff),-1*x_cent,-1*y_cent)
endfor

return,ccube

end


