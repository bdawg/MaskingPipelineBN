; split a frame into quadrants with the standard
; pharo ordering (somewhat arcane)
; $Id: pharosplit.pro,v 1.1.1.1 2005/12/19 05:15:05 mireland Exp $

function pharosplit, frame
   sz = size(frame)

   xs = sz[1]/2
   ys = sz[2]/2
 
   x0 = 0
   x1 = xs-1
   x2 = xs
   x3 = 2*xs-1

   y0 = 0
   y1 = ys-1
   y2 = ys
   y3 = 2*ys-1

   quadrants = fltarr(xs,ys,4)
   quadrants[*,*,1] = frame[x0:x1,y0:y1]
   quadrants[*,*,0] = frame[x2:x3,y0:y1]
   quadrants[*,*,2] = frame[x0:x1,y2:y3]
   quadrants[*,*,3] = frame[x2:x3,y2:y3]

   return, quadrants
end

$Log: pharosplit.pro,v $
Revision 1.1.1.1  2005/12/19 05:15:05  mireland
Initial release of masking code.

Revision 1.1  2003/12/31 02:20:02  jpl
continued development, including flatfields and bad pixel treatment
for pharoreadfits

