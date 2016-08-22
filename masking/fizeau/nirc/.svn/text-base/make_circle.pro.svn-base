pro make_circle,im,x,y,r,sz,pvect=pvect,ellipse=ellipse,double=double
; Program makes a circle at given location in image.
; Written by John sometime
; PGT - added option to pass out pixel vector PVECT
; ACC - added option to make an ellipse instead of a circle
; ellipse= [semimajor axis,semiminor axis]

if (keyword_set(im) eq 0) then begin
  print,'pro make_circle,im,x,y,r,sz'
  return
endif

if (keyword_set(sz) eq 0) then sz=1.0
info=size(im)

if (keyword_set(ellipse) eq 0) then ellipse=[1.,1.]
if (n_elements(ellipse) ne 2) then begin
    print,'to make an ellipse you need a semimajor and semiminor axis'
    stop
endif

a=(ellipse[0])^2
b=(ellipse[1])^2

xxx=dindgen(info(1))
proj=replicate(1,info(1))
xx=xxx#proj
yy=proj#xxx
xx=xx-x
yy=yy-y

IMage=sqrt(1/a*xx*xx+1/b*yy*yy)
pvect=where(image le r)
image(pvect)=1.0*sz
image(where(image gt r))=0.0
if keyword_set(double) eq 0 then image = float(image)
IM=IM+image

;for i=0,info(1)-1 do begin
;for j=0,info(2)-1 do begin
;   r_d=sqrt((float(i)-float(x))^2+(float(j)-float(y))^2)
;   if (r_d le r) then im(i,j)=im(i,j)+size
;endfor
;endfor
end

