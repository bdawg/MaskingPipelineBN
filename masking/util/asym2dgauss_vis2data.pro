;2003Sep28 JDM
;2005Nov20 JDM
;2005Nov23 JDM

function asym2dgauss_vis2data,params,vis2data=vis2data,file=file,$
  scale=scale,dim=dim,im=im,x=x,y=y,aa=aa,tt=tt,squarex=squarex,$
  squarey=squarey

if keyword_set(dim) eq 0 then dim=64
if keyword_set(scale) eq 0 then scale = 1.0 ;mas per pix

; For use with oidata library by acting on vis2data from
; extract_vis2data
; alternatively, can take a oidata file directly!
; INPUTS:
; Model 
; Params:
;params=[point_flux,gauss FWHM major, gauss FWHM MINOR, Major Angle,skew,theta0-deg]



if (keyword_set(vis2data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(vis2data) eq 0) then begin
  extract_vis2data,file=file,vis2data
endif

;else assume we got vis2data

model_vis2data=vis2data
im=fltarr(dim,dim)
modparams=params(1:6)
modparams(0)=modparams(0)*(1./(params(0)+params(1)))

makeasym2dgauss,im,modparams,scale=scale,x=x,y=y,aa=aa,tt=tt,$
   squarex=squarex,squarey=squarey
im=im*(params(0)+params(1))
;stop
extract_data,im,visib,phases,scale=scale,u=vis2data.u,v=vis2data.v,/cubic,$
  /nohan

;uniform_disk,sqrt(vis2data.u*vis2data.u+vis2data.v*vis2data.v),[mas2rad(params(1)),params(0)],ampv
;
; POINT SOURCE built-in to makring

visib_real =      visib*cos(phases*!pi/180.)
visib_imag =      visib*sin(phases*!pi/180.)

model_vis2data.vis2data=visib_real^2 + visib_imag^2
model_vis2data.vis2err=0.
;stop
return,model_vis2data
end




