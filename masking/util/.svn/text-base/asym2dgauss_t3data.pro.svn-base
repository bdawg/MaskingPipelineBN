;2003Sep28 JDM
;2005Nov20 JDM
;2005Nov23 jdm


function asym2dgauss_t3data,params,t3data=t3data,file=file,scale=scale,dim=dim,$
 im=im,x=x,y=y,aa=aa,tt=tt,squarex=squarex,squarey=squarey


; For use with oidata library by acting on t3data from
; extract_t3data
; alternatively, can take a oidata file directly!

; INPUTS:
; Model
; Params:
;params=[point_flux,gauss FWHM major, gauss FWHM MINOR, Major Angle,skew,theta0-deg]

if (keyword_set(t3data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(t3data) eq 0) then begin
  extract_t3data,file=file,t3data
endif

;else assume we got  t3data

model_t3data=t3data

im=fltarr(dim,dim)
modparams=params(1:6)
modparams(0)=modparams(0)*(1./(params(0)+params(1)))

makeasym2dgauss,im,modparams,scale =scale,x=x,y=y,aa=aa,tt=tt,$
  squarex=squarex,squarey=squarey  ; makes total in image 1.0
im=im*(params(0)+params(1))
;stop
extract_data,im,visib1,phases1,scale=scale,u=t3data.u1,v=t3data.v1,/cubic,$
  /nohan
;uniform_disk,sqrt(t3data.u1*t3data.u1+t3data.v1*t3data.v1),$
;   [mas2rad(params(1)),params(0)],ampv1
; POINT SOURCE built-in to makring

visib1_real =      visib1*cos(phases1*!pi/180.)
visib1_imag =      visib1*sin(phases1*!pi/180.)

extract_data,im,visib2,phases2,scale=scale,u=t3data.u2,v=t3data.v2,/cubic,$
  /nohan
;uniform_disk,sqrt(t3data.u2*t3data.u2+t3data.v2*t3data.v2),$
;   [mas2rad(params(1)),params(0)],ampv2
visib2_real =      visib2*cos(phases2*!pi/180.)
visib2_imag =      visib2*sin(phases2*!pi/180.)

extract_data,im,visib3,phases3,scale=scale,u=t3data.u3,v=t3data.v3,/cubic,$
  /nohan
;uniform_disk,sqrt(t3data.u3*t3data.u3+t3data.v3*t3data.v3),$
;   [mas2rad(params(1)),params(0)],ampv3
visib3_real =      visib3*cos(phases3*!pi/180.)
visib3_imag =      visib3*sin(phases3*!pi/180.)

ri2at,visib1_real,visib1_imag,visib1,phases1
ri2at,visib2_real,visib2_imag,visib2,phases2
ri2at,visib3_real,visib3_imag,visib3,phases3


model_t3data.t3amp=visib1*visib2*visib3
model_t3data.t3amperr=0.
model_t3data.t3phi=mod360(phases1+phases2+phases3)
model_t3data.t3phierr=0.
;stop
return,model_t3data
end




