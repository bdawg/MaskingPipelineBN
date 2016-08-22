;2003Sep28 JDM


function binary_t3data,binary_params,t3data=t3data,file=file

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

if (keyword_set(t3data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(t3data) eq 0) then begin
  extract_t3data,file=file,t3data
endif

;else assume we got  t3data

model_t3data=t3data

binary_disks,t3data.u1,t3data.v1,binary_params,visib1,phases1
binary_disks,t3data.u2,t3data.v2,binary_params,visib2,phases2
binary_disks,t3data.u3,t3data.v3,binary_params,visib3,phases3

model_t3data.t3amp=visib1*visib2*visib3
model_t3data.t3amperr=0.
model_t3data.t3phi=mod360(phases1+phases2+phases3)
model_t3data.t3phierr=0.

return,model_t3data
end




