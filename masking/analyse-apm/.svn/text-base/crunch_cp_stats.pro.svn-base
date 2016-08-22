;;
;; Calculate mean, error, covariance for closure phase and bispectrum
;; amplitude
;;
;;  NOTE: Returns CPs in DEGREES
;;
;;  Author: David Bernat, dwb29@cornell.edu

pro crunch_cp_stats, bs_arr, stats

  sz = size( bs_arr, /dimension )
  n_arr = sz[0]
  n_bispect = sz[1]

;;  Now calculate the mean and stddev of these bispectrum
;;  This routine calculates the average and standard deviation of
;;  closure phases directly, rather than the RMS of the bispectrum

;; To avoid the discontinuity at 180 degrees (ie 180 + 1 = -179),
;; rotate the bispectrum so their mean is 0.

  bs_mean = total( bs_arr, 1 ) / n_arr
  bs_rot = bs_arr
  For i = 0, n_bispect - 1 DO bs_rot[*,i] *= conj( bs_mean[i] ) / abs( bs_mean[i] )

  cp_mean = ( atan( bs_mean, /phase ) + total( atan( bs_rot, /phase ), 1 ) / n_arr ) * !radeg
  cp_arr_meansub = atan( bs_arr, /phase ) * !radeg - Transpose( temporary( rebin( cp_mean, n_bispect, n_arr, /sample ) ) )

  cp_cov = ( Transpose( cp_arr_meansub ) # ( cp_arr_meansub ) ) / n_arr

  cp_err = dblarr( n_bispect )
  for i = 0, n_bispect - 1 DO cp_err[i] = Sqrt( cp_cov[i,i] )

  amp2_mean = ABS( total( bs_arr, 1 ) / n_arr )^2
  amp2_arr_meansub = ABS( bs_arr )^2 - Transpose( temporary( rebin( amp2_mean, n_bispect, n_arr, /sample ) ) )
  
  amp2_cov = Transpose( amp2_arr_meansub ) # amp2_arr_meansub / n_arr
  amp2_err = dblarr( n_bispect )
  for i = 0, n_bispect - 1 DO amp2_err[i] = Sqrt( amp2_cov[i,i] )

  stats = { cp_mean: cp_mean, cp_err: cp_err, cp_cov: cp_cov, amp2_mean: amp2_mean, amp2_err: amp2_err, amp2_cov: amp2_cov }

END
