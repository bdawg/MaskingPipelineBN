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


function calc_loglike, cps, datalikes
  res = fltarr( n_elements( cps ) )
  For i = 0, n_elements( cps ) - 1 DO res[i] = alog( (interpol( (*datalikes[i]).l_vec, (*datalikes[i]).m_vec, cps[i] ) > 1D-200 ) < 1.0 )
  return, Sqrt( -res )
END

; Function used for MPFIT in my Binary Grid Script
function mp_binary_func_fast, params, dp, t3data=t3data, EVecs=Evecs, t3template=t3template, datalikes=datalikes, use_angle_diff=use_angle_diff
t3model = binary_t3data_fast(params, t3data=t3template)

;; t3model returns multivariate Closure Phases, Evecs maps them to
;; their linearly independent sets.  (Evecs are the eigenvectors of
;; the covariance matrix)

If( Keyword_Set( Evecs ) ) THEN $
   residuals = [( t3data.t3phi - Evecs # t3model.t3phi) / t3data.t3phierr ] $
ELSE If( Keyword_Set( datalikes ) ) THEN $
   If( use_angle_diff ) $
   THEN residuals = calc_loglike( angle_diff( t3model.t3phi, t3data.t3phi ), datalikes ) $
   ELSE residuals = calc_loglike( t3model.t3phi, datalikes ) $
ELSE IF( Keyword_Set( use_angle_diff ) ) THEN $
   residuals = angle_diff(t3data.t3phi, t3model.t3phi)/t3data.t3phierr $
ELSE $
   residuals = (t3data.t3phi - t3model.t3phi)/t3data.t3phierr

return, residuals
end

; Calculate Chi2  (non-reduced, now)
function binary_chi2_many, params, dp, t3model=t3model, t3data=t3data, datalikes=datalikes, use_angle_diff=use_angle_diff

  sz = size( t3model, /dimension )
  n_wide = sz[1]
  chi2 = dblarr( n_wide )

  If( keyword_set( datalikes ) ) THEN BEGIN
     If( keyword_set( use_angle_diff ) ) THEN BEGIN      
        
        For i = 0, n_wide - 1 DO chi2[i] = total( (calc_loglike( angle_diff( (t3model.t3phi)[*,i], t3data.t3phi ), datalikes ) )^2 )

     ENDIF ELSE BEGIN

        For i = 0, n_wide - 1 DO chi2[i] = total( (calc_loglike( (t3model.t3phi)[*,i], datalikes ))^2 )

     ENDELSE
  ENDIF ELSE BEGIN

     ;; Covariances / Gaussian Errors
     If( keyword_set( use_angle_diff ) ) THEN BEGIN
        For i = 0, n_wide - 1 DO chi2[i] = total( (angle_diff(t3data.t3phi, (t3model.t3phi)[*,i])/t3data.t3phierr)^2 )
     ENDIF ELSE BEGIN

        For i = 0, n_wide - 1 DO chi2[i] = total( ((t3data.t3phi - (t3model.t3phi)[*,i] )/t3data.t3phierr )^2 )
     ENDELSE
  ENDELSE

  return, chi2
end
