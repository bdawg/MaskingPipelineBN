;PGT Nov02 - 	added passout of origin of image1(0,0) in x,y 
;          -    Fixed bug (did not work for odd n)
;          
;   if x,y passed in then the program grabs a region centered
;          on this. x,y are OVERWRITTEN however with the location
;          of the origin of the new array in coords of the input image. 

function grabnxn,image,n,x=x,y=y,wrap=wrap

if (not keyword_set(x) or not keyword_set(y)) then begin
  if (min(image) lt max(image)) then locate_peak,image,x,y else begin
		x=(size(image))(1)/2.
		y=(size(image))(2)/2.
		endelse
endif

dif=n/2

if (keyword_set(wrap)) then begin
; Case to preserve size by wrapping array if it is near edge
  image1 = shift(image,dif-x,dif-y)
  image1=image1[0:n-1,0:n-1]

endif else begin
; Two cases for even or odd number of pixels in map

  if( dif*2 ne n) then begin	; Case1: ODD n
    xmin=x-dif > 0
    xmax=x+dif < ( (size(image))(1)-1 )
    ymin=y-dif > 0
    ymax=y+dif < ( (size(image))(2)-1 )
    image1=image(xmin:xmax,ymin:ymax)
  endif else begin		; Case2: EVEN n
    xmin=x-dif   > 0
    xmax=x+dif-1 < ( (size(image))(1)-1 )
    ymin=y-dif   > 0
    ymax=y+dif-1 < ( (size(image))(2)-1 )
    image1=image(xmin:xmax,ymin:ymax)
  endelse
  x=xmin
  y=ymin
endelse

return,image1
end


