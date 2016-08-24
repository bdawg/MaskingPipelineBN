function fit_mf1_fitfunc, x, y, p
  ; Make a proposed mf image
  
common commonblock
restore,inputMF
if flip eq 1 then xy_coords[*,0] = -xy_coords[*,0]

xscale = p[0]
yscale = p[1]
msk_rot = p[2]
;; scale = p[0]
;; msk_rot = p[1]

rot1=[[cos(msk_rot),sin(msk_rot)],[-sin(msk_rot),cos(msk_rot)]]
;; xy_coords = scale*rot1##xy_coords

;; rot1[*,0] = xscale*rot1[*,0]
;; rot1[*,1] = yscale*rot1[*,1]
xy_coords = rot1##xy_coords
xy_coords[*,0] = xscale*xy_coords[*,0]
xy_coords[*,1] = yscale*xy_coords[*,1]
;;xy_coords = rot1##xy_coords

r=makeSplodgeR
mfMask = dblarr(chipsz,chipsz)

for bl=0,n_baselines-1 do begin
     thisholepair=bl2h_ix[*,bl]
     Xposn =-(xy_coords[bl2h_ix[0,bl],0] - xy_coords[bl2h_ix[1,bl],0])/filter[0] $
           *rad_pixel*chipsz+chipsz/2
     Yposn =-(xy_coords[bl2h_ix[0,bl],1] - xy_coords[bl2h_ix[1,bl],1])/filter[0] $
           *rad_pixel*chipsz+chipsz/2

;     gaussparams = [ampl, sampledisk_r, 0, Xposns[bl], Yposns[bl] ]
;     currentSplodge = gauss_circ(chipsz, gaussparams)

     d = sqrt( (xposn-chipsz/2)^2 + (yposn-chipsz/2)^2) / (chipsz/2)
     wt = 1. - longBLweighting*d^2

     currentSplodge = wt*maskAmpl*exp(-( (x-Xposn)^2/(2.*r^2) + $
                               (y-Yposn)^2/(2.*r^2) ))
     mfMask = mfMask + currentSplodge
endfor

if showimage then begin
   ;image_cont,logps-mfmask,/noc,/asp
   sz = (size(logps))[1]
   tvscl,rebin((logps-mfmask)>(-3),sz*3,sz*3)
endif

return, mfMask
end
