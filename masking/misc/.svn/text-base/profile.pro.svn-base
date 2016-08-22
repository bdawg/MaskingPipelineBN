pro profile,array


; JDM 1/20/97

; profile,image
;
; This routine will display an image, then let you choose two points 
; and then plot the profile beneath it in a separate window.
;
; Just click out of bounds to start over.
;
x1=0
y1=0
x2=0
y2=0

info=size(array)
window,1,xsize=500,ysize=300
window,0,xsize=500,ysize=500
image_cont,array,/asp,/nocont
startagain:
wset,0
print,'Click on two point with left mouse button:'
image_cont,array,/asp,/nocont
plots,[x1,x2],[y1,y2]
cursor,x1,y1 & wait,.5

if (x1 lt 0 or  x1 ge info(1) or y1 lt 0 or y1 ge info(2) ) then return
cursor,x2,y2
if (x2 lt 0 or  x2 ge info(1) or y2 lt 0 or y2 ge info(2) ) then return

plots,[x1,x2],[y1,y2]

wset,1
xval=findgen(100)/99.
yval=findgen(100)/99.
xval=xval*(x2-x1)+x1
yval=yval*(y2-y1)+y1
;stop
plot,prof(array,xval,yval)
wait,.5
goto,startagain
end



