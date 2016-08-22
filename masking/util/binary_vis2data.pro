;2003Sep28 JDM


function binary_vis2data,binary_params,vis2data=vis2data,file=file

; For use with oidata library by acting on vis2data from
; extract_vis2data
; alternatively, can take a oidata file directly!

; INPUTS:
; Model of a binary star, each with UD sizes and a ratio.
; Params:
;   params(0) = Separation (mas)
;   params(1) = Position Angle (degs,
;               E of N, pointing from primary -> secondary)
;   params(2) = Brightness Ratio of Primary over Secondary
;   params(3) = UD size of primary (mas)
;   params(4) = UD size of secondary (mas)

if (keyword_set(vis2data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(vis2data) eq 0) then begin
  extract_vis2data,file=file,vis2data
endif

;else assume we got vis2data

model_vis2data=vis2data

binary_disks,vis2data.u,vis2data.v,binary_params,visib,phases

model_vis2data.vis2data=visib^2
model_vis2data.vis2err=0.

return,model_vis2data
end




