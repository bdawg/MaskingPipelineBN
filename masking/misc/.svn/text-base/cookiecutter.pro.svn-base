pro cookiecutter,im,x,y,r,sz,pvect=pvect
; Program makes a set of circles at specified locations in image.
; Takes input arrays x and y and makes circles of value sz 
; (default 1) in the image im

if (keyword_set(im) eq 0) then begin
  print,'pro cookiecutter,im,x,y,r,sz'
  return
endif

if (keyword_set(sz) eq 0) then sz=1.0
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

