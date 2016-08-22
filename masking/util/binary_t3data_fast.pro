;; ----------------------------------------------------
;;  Binary_t3data_fast.pro and binary_chi2_fast.pro
;;   Both these files have been written to obtain model t3data or
;;   calculate chi2 as quickly. This is helpful for monte carlo
;;   simulations that call these functions many many times.
;;
;;  NOTE: Both *_fast files remove parameters 3 and 4 (the angular
;;  size of the stars and the stars are appropriately treated as point
;;  sources.  IE, parms = [ 100., 90., 150. ]
;;  
;;  Author: David Bernat, dwb29@cornell.edu
;;




function binary_t3data_fast, binary_params, t3data=t3data

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

model_t3data=t3data

U = [ t3data.u1, t3data.u2, t3data.u3 ]
V = [ t3data.v1, t3data.v2, t3data.v3 ]

binary_fast, U, V, binary_params, phases

n = n_elements( t3data.u1 )
phases = reform( phases, n, 3 )

model_t3data.t3phi=mod360( total( phases, 2 ) )
model_t3data.t3phierr=0.

return,model_t3data
end

