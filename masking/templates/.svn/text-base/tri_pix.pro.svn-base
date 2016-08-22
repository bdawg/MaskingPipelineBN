; Function to return a 3 by ncpx array which indexes all closing
; triplets of pixels. The 3 splodges are all assumed to lie at the 
; center (npix/2) of an array of dimension npix. When they are used
; a simple shift of the splodge can be implemented by adding or subtracting
; from the pixel locations.

function tri_pix, array_size, hole_radius  

d=fltarr(array_size,array_size)
cookiecutter,d,array_size/2,array_size/2,hole_radius

pvct=where(d gt 0)
npx=n_elements(pvct)

; Now build array of valid CLP triangles for this single splodge
;px1=0

for px1=0,npx-1 do begin
  thispix1=array_coords(pvct[px1],d)
  xcor_12=d+shift(d,array_size/2-thispix1[0],array_size/2-thispix1[1])
  valid_b12_vct=where(xcor_12 gt 1)
  thisntriv=n_elements(valid_b12_vct)
  thistrivct=lonarr(3,thisntriv)
  ; Now loop over the number of valid px2 given starting px1, accumulate triangles
  for px2=0,thisntriv-1 do begin
    ; now find the right xy_coord to close the triangle given bl1 and bl2 ...
      thispix2=array_coords(valid_b12_vct[px2],d)
      thispix3=array_size*1.5-(thispix1 + thispix2)
      bl3_px=thispix3[1]*array_size+(thispix3[0])
      thistrivct[*,px2]=[pvct[px1],valid_b12_vct(px2),bl3_px]
  endfor
  if(px1 eq 0) then closing_tri_pix=thistrivct else closing_tri_pix=[[closing_tri_pix],[thistrivct]]
endfor

return,closing_tri_pix

end
