;This is a program to azimuthally average an array

function azi_ave,  a,  av2d = av2d

sz =  (size(a))[1]
d =  shift(dist(sz), sz/2,  sz/2)
nelt =  fltarr(sz)
for i =  long(0),n_elements(d)-1 do nelt[d[i]]++
azi_ave =  fltarr(sz)
for i =  long(0),n_elements(d)-1 do azi_ave[d[i]] += a[i]
azi_ave =  azi_ave/(nelt > 1)
if (arg_present(av2d)) then begin
 av2d =  a
 av2d[*] = azi_ave[d[*]]
endif

return,  azi_ave

end
