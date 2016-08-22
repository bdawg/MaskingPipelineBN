;This function returns a pupil with specified properties
function make_pupil, dim, diam, focus=focus, a1=a1

if (keyword_set(focus) eq 0) then focus = 0
if (keyword_set(a1) eq 0) then a1 = 0 
pupil = complexarr(dim,dim)
mini = long(dim/2 - diam/2)
maxi = dim/2 + diam/2
rad = diam/2.0
for i = mini,maxi do for j = mini,maxi do begin
 x = dim/2.0 - i
 y = dim/2.0 - j
 dist2 = (dim/2-i)^2 + (dim/2-j)^2
 if (dist2 lt rad^2)  then $
  pupil[i,j] = exp(complex(0,x^2/rad^2*(focus -a1) + y^2/rad^2*(focus+a1)))
endfor

return, pupil
end
