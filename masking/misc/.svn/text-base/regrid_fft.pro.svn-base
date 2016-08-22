;This function regrids an array by an fft, zero-padding and an
;inverse-fft. The array is assumed to be periodic unless xrev or yrev
;are set.

function regrid_fft, inarr,  dimx,  dimy, xrev = xrev,  yrev = yrev
a = inarr
inx =  (size(a))[1]
iny =  (size(a))[2]
if (dimx lt inx or dimy lt iny) then stop
if (keyword_set(xrev)) then begin
 b =  fltarr(2*inx,  iny)
 b[0:inx-1, *] =  a
 b[inx:2*inx-1, *] = reverse(a, 1)
 a =  b
 inx *= 2 
 dimx *= 2
endif
if (keyword_set(yrev)) then begin
 b =  fltarr(inx,  2*iny)
 b[*, 0:iny-1] =  a
 b[*, iny:2*iny-1] = reverse(a, 2)
 a =  b
 iny *= 2 
 dimy *= 2
endif

outft =  complexarr(dimx, dimy)
outft[0:inx-1, 0:iny-1] =  shift(fft(a, -1), inx/2, iny/2)
outft =  shift(outft, -inx/2, -iny/2)
out = float(fft(outft, 1))

if (keyword_set(xrev)) then dimx /= 2
if (keyword_set(yrev)) then dimy /= 2

out =  out[0:dimx-1, 0:dimy-1]

return,  out

end
