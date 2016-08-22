;;----------------------------------------
;;
;;  The _many files are modified versions of the _fast files to take
;;   advantage of IDL's ultrafast processing of arrays.  These
;;   files are adapted from some files originally from John Monnier.  See
;;   the _fast files for more notes.
;;
;; NOTE: parms 3 and 4 (angular size of stars) have been removed for
;;  speed considerations.  See _fast files for more details.
;;
;; NOTE: Returns an array of Model CPs only, and not an array of
;;  t3data structures filled with Model CPS
;; Author: David Bernat, dwb29@cornell.edu
;;
;; ----------------------------------------


function binary_t3data_many, binary_params, t3data=t3data, Return_CPs_Only=Return_CPs_Only

  If( NOT Keyword_Set( Return_CPs_Only ) ) THEN Return_CPs_Only = 0

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

model_t3data=t3data

sz = size( binary_params, /dimensions ) 
n_wide = sz[1]

sz = size( t3data.u1, /dimensions )
n_u = sz[0]

;; Parallelize
U = rebin( [ t3data.u1, t3data.u2, t3data.u3 ], 3*n_u, n_wide, /sample )
V = rebin( [ t3data.v1, t3data.v2, t3data.v3 ], 3*n_u, n_wide, /sample )

binary_many, U, V, binary_params, phases

phases = reform( phases, n_u, 3, n_wide )

If( Return_CPs_Only EQ 0 ) THEN BEGIN
   For i = 0, n_wide - 1 DO $
      If( i EQ 0 ) THEN model_t3data = [ t3data ] ELSE model_t3data = [ [ model_t3data ], [t3data ] ] 
   model_t3data.t3phi=mod360( total( phases, 2 ) )
   model_t3data.t3phierr=fltarr( n_u, n_wide)
   Return, model_t3data
ENDIF ELSE return, mod360( total( phases, 2 ) )

End

