;;----------------------------------------
;;
;;  The _paral files are modified versions of the _fast files to take
;;   advantage of IDL's ultrafast processing of arrays.  These
;;   files are adapted from some files originally from John Monnier.  See
;;   the _fast files for more notes.
;;
;; NOTE: parms 3 and 4 (angular size of stars) have been removed for
;;  speed considerations.  See _fast files for more details.
;;
;; Author: David Bernat, dwb29@cornell.edu
;;
;; ----------------------------------------


function binary_t3data_paral, binary_params, t3data=t3data

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

binary_paral, U, V, binary_params, phases

s = size( t3data.u1, /dimensions)
n_u = s[0]
n_wide = s[1]

phases = reform( phases, n_u, 3, n_wide )
model_t3data.t3phi=mod360( total( phases, 2 ) )
model_t3data.t3phierr=fltarr( n_u, n_wide)

return,model_t3data
end

