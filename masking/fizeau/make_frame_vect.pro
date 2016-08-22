;######## RELEASE CLIP HERE ######## 
;This procedure minimizes effort in running the cube-making software
;for nirc2 or pharo.

function make_frame_vect,  fstart,  nfs,  tsize = tsize

if (n_elements(fstart) ne n_elements(nfs)) then begin
  print,  'Error! fstart and nfs must have the same number of elements...'
  stop
endif
frames = -1
tsize = -1
for i = 0, n_elements(fstart)-1 do begin
 frames =  [frames, fstart[i]+indgen(abs(nfs[i]))]
 tsize =  [tsize, replicate(sign(nfs[i])*(i+1.)/10., abs(nfs[i]))]
endfor
frames =  frames[1:*]
tsize =  tsize[1:*]

return,  frames

end
