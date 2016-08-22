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
;; NOTE: Now returned UNREDUCED Chi2
;;
;; Author: David Bernat, dwb29@cornell.edu
;;
;; ----------------------------------------

; Function used for MPFIT in my Binary Grid Script

function mp_binary_func_paral, params, dp, t3data=t3data
t3model = binary_t3data_paral(params, t3data=t3data)
residuals = [mod360(t3data.t3phi - t3model.t3phi)/t3data.t3phierr]
return, residuals

end


;
; Calculate Chi2  (non-reduced, now)
;

function binary_chi2_paral, params, dp, t3model=t3model, t3data=t3data
call_it_chi = mod360(t3data.t3phi - t3model.t3phi)/t3data.t3phierr
chi2 = total(call_it_chi^2, 1 )
return, chi2

end
