;;----------------------------------------
;;
;;  The _many files are modified versions of the _fast files to take
;;   advantage of IDL's ultrafast processing of arrays.  These
;;   files are adapted from some files originally from John Monnier.  See
;;   the _fast files for more notes.
;;
;; Author: David Bernat, dwb29@cornell.edu
;;
;; ----------------------------------------


pro binary_many, U, V, params, phases

delta_dec = mas2rad( params[0, *] * Cos( params[1, *] / !radeg ) )
delta_ra  = mas2rad( params[0, *] * Sin( params[1, *] / !radeg ) )
;recall +ra points at pa 90

s = size( U, /dimensions)
n_u = s[0]
n_wide = s[1]

delta_dec = rebin( delta_dec, n_u, n_wide, /sample ) 
delta_ra  = rebin( delta_ra,  n_u, n_wide, /sample ) 

; Using boden notation (I hope!) from michelson book
arg = 2.0 * !pi * ( u * delta_ra + v * delta_dec )

; No reason to normalize flux (see binary_disks) since we only care
; about phase

primary_factor = rebin( params[2,*], n_u, n_wide, /sample )
complex_vis = primary_factor + complex( Cos(arg), Sin( arg) )
phases = atan( complex_vis, /phase ) * !radeg

end
