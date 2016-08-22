function ave_azi,image

ratio=1e8/max(image)
image=image*ratio
; New PGT version: no need for square image...
info=size(image)
nx=info(1)
ny=info(2)

distance=float(dist(nx,ny))

results=histogram(distance,r=r)
aveazi=float(results)

for i=0,n_elements(results)-1 do begin
         aveazi(i)=total(image(r(r(i):r(i+1)-1)))/float(results(i))
endfor
image=image/ratio
return,aveazi/ratio
end
