;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;      cookiecutter_conica
;
; PURPOSE:
;      Overwrites fixed-radius circular regions at specified locations in a 2D array 
;      with a specified value. This is useful for constructing 'cookiecutters' 
;      i.e. 2D arrays with circular regions of value 0 and value 1 everywhere else.      
;
; CALLING SEQUENCE:
;      cookiecutter_conica,im,x,y,r,sz
;
; INPUTS:
;      im - input 2D array
;      x - vector containing x-coordinates of circles to be cut out
;      y - vector containing y-coordinates of circles to be cut out
;      r - scalar value specifying radius of circles to be cut out
;      sz - value to be used for the cut out circle regions
;
; OUTPUT:
;      im - output 2D array
;
; HISTORY:
;      version 1.0  TME 28May09  Adapted from cookiecutter.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro cookiecutter_conica,im,x,y,r,sz


if (keyword_set(im) eq 0) then begin
  print,'pro cookiecutter,im,x,y,r,sz'
  return
endif

info=size(im)
s=size(x)

; Check to see if we have been given a vector or scalar
; array of holes to cut.
if(s(0) lt 1) then begin
    xx=[x]
    yy=[y] 
    ncirc=1
endif else begin
    ncirc=s(1)
    xx=x
    yy=y
endelse

for c=0,ncirc-1 do begin
   for i=fix(max([0.0,xx(c)-r-1])),fix(min([info(1)-1,xx(c)+r+1])) do begin
   for j=fix(max([0.0,yy(c)-r-1])),fix(min([info(2)-1,yy(c)+r+1])) do begin
      r_d=sqrt((float(i)-float(xx(c)))^2+(float(j)-float(yy(c)))^2)
      if (r_d le r) then im(i,j)=sz
   endfor
   endfor
endfor

end

