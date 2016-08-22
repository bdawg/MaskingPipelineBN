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
;;  NOTE: This returns UNREDUCED Chi2
;;  
;;  Author: David Bernat, dwb29@cornell.edu
;;

; Function used for MPFIT in my Binary Grid Script
function mp_binary_func_fast, params, dp, t3data=t3data
t3model = binary_t3data_fast(params, t3data=t3data)
residuals = [mod360(t3data.t3phi - t3model.t3phi)/t3data.t3phierr]
return, residuals
end

;
; Calculate Chi2  (non-reduced, now)
;

function binary_chi2_fast, params, dp, t3model=t3model, t3data=t3data
call_it_chi = mod360(t3data.t3phi - t3model.t3phi)/t3data.t3phierr
chi2 = total(call_it_chi^2)
return, chi2

end
