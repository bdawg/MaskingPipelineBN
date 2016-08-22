;This procedure fits a function of the form:
;V^2 = (a0 + a1*r + a2/r)*binary_vis2
;... to the input vis2data.
function v2_2nd_order,  X,  Y,  P

return,  p[0] + p[1]*X^2 + p[2]*Y^2 + p[3]*X*Y
end

function binary_vis2atm,binary_params,vis2data=vis2data, atm = atm,  assym = assym

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

if (keyword_set(vis2data) eq 0 ) then begin
 print,'Must input some data!'
return,-1
endif


model_vis2data=vis2data

binary_disks,vis2data.u,vis2data.v,binary_params,visib,phases

if (keyword_set(assym)) then begin
  u =  vis2data.ucoord
  v =  vis2data.vcoord
  p =  mpfit2dfun('v2_2nd_order', u,  v,  vis2data.vis2data/visib^2,  $
                  vis2data.vis2err/visib^2, [1.0, 0, 0, 0],  /quiet,  yfit = atm)
endif else begin
  r =  sqrt(vis2data.ucoord^2+ vis2data.vcoord^2)
  expr =  'p[0] + p[1]*X'
  p =  mpfitexpr(expr, r,  vis2data.vis2data/visib^2,  vis2data.vis2err/visib^2, [1.0, 0],  /quiet,  yfit = atm)
endelse

model_vis2data.vis2data = visib^2*atm
model_vis2data.vis2err=0.

return,model_vis2data
end

