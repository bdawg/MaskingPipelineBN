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
;; Adapted from JDM's binary_disks program.

pro binary_fast, U, V, params, phases

delta_dec = mas2rad( params(0) * Cos( params[1] / !radeg ) )
delta_ra  = mas2rad( params(0) * Sin( params[1] / !radeg ) )
;recall +ra points at pa 90

; Using boden notation (I hope!) from michelson book
arg = 2.0 * !pi * ( u * delta_ra + v * delta_dec )

; No reason to normalize flux (see binary_disks) since we only case
; about phase
complex_vis = 1.0 + params(2) * complex( Cos(arg), Sin( arg) )
phases = atan( complex_vis, /phase ) * !radeg

end
