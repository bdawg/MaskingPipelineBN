;; ----------------------------------------------
;;  This code uses a bootstrapping method to calculate the error
;;   and covariances between closure phases from bispectrum.
;;   Bootstrapping is a more accurate method at low signals to noise
;;   (errors tend to increase by ~ 50% using bootstrapping over a
;;   gaussian approximation).  
;;
;; This structure of this code is written to make use of IDLs
;; ultrafast processing of arrays and its distain for nested loops.
;; The rn array is very large (in fact, capped so as to not overload
;; the RAM).  Modifications to this code should be done with care (but
;; improvements are of course invited.)
;;
;;  Author: David Bernat, dwb29@cornell.edu
;;


pro bootstrap_cp_err, bispect_targ, bispect_cals, cp_cov, cp_err, cp_mean, v2_cov, v2_err, v2_mean, bs_arr, nSim=nSim

  If( NOT keyword_set( nSim ) ) THEN nSim = 3000
 
  sz = size( bispect_targ, /dimensions )
  n_targ    = sz[0] ; Number of data points
  n_bispect = sz[1] ; Number of bispectrum per exposure (usually 84)

  sz = size( bispect_cals, /dimensions )
  n_cals    = sz[0] ; Number of data points
  n_bispect2 = sz[1] ; Number of bispectrum per exposure (usually 84)

  ; I want this line to crash since this runs in a cron
  ; If( n_bispect NE n_bispect2 ) THEN stop

  ; This is a semi-arbitrary number that I pulled from my own anecdotal
  ; testing.  This is the number of elements in the operational array.
  ; Too large and it crushes the RAM. (1 float = 4 bytes)
  limit = 25.e6

  n_sim_per = Ceil( limit / (n_targ + n_cals ) / n_bispect ) < nSim
  n_bunches = Ceil( 1.0*nSim / n_sim_per )

  nSim = n_bunches * n_sim_per
  seeds_targ = 2L + lindgen( n_bunches )
  seeds_cals = 200L - lindgen( n_bunches )

  cp_arr = fltarr( nSim, n_bispect )
  bs_arr = complexarr( nSim, n_bispect )
  bs_targ_arr = complexarr( n_sim_per, n_bispect )
  bs_cals_arr = complexarr( n_sim_per, n_bispect )
  print, "Starting Bootstrapping Simulation..."
  For n = 0, n_bunches - 1 DO BEGIN
     rn_targ = floor( randomu( seeds_targ[n], n_sim_per * n_targ ) * n_targ )
     rn_cals = floor( randomu( seeds_cals[n], n_sim_per * n_cals ) * n_cals )

     ; This used to be done on a per-bispectrum level -- otherwise
     ; the array sizes may be large enough to crash the RAM.  But
     ;  now the bispectrum arrays are on the order of a few hundred, not a
     ;  few ten thousand.

     ;; Now, randomly select single bispect data to mock up a new
     ;; measurement of size n_data
     ;; NOTE: This is an extremely memory intensive process.  Breaking
     ;; this line up into its composite parts requires additional
     ;; variables to be created and kept around...ie, memory leaks.
     ;; NOTE: The RHS is an array, the LHS appears to be a single
     ;; float variable (notice the indices).  However, IDL is smart
     ;; enough to use [0,i] as the starting point for the RHS'
     ;; array, effectively filling in [*,i] as the dimensions match.
     ;; You can thank DFanning for that one.

        bs_targ_arr = total( reform( temporary( bispect_targ[rn_targ,*] ), n_sim_per, n_targ, n_bispect, /overwrite ), 2 )
        bs_cals_arr = total( reform( ( temporary( bispect_cals[rn_cals,*]) ), n_sim_per, n_cals, n_bispect, /overwrite ), 2 )
        bs_arr[(n*n_sim_per):((n+1)*n_sim_per-1),*] = bs_targ_arr / bs_cals_arr
        cp_arr[(n*n_sim_per):((n+1)*n_sim_per-1),*] = atan( bs_targ_arr / bs_cals_arr, /phase )

  print, "Finished Bunch ", n+1, " of ", n_bunches, " bunches.  "
  print, "Finished ", (n+1)*n_sim_per, " bootstraps." 
  ENDFOR

  ;;  Now calculate the stddev of these bispectrum
  ;;   Notice that I use the cp_mean of the input bispectrum, and not
  ;;   the mean of the bootstrapped data.
  cp_mean = atan( total( bispect_targ, 1 ) / total( bispect_cals, 1 ), /phase )
  cp_arr_meansub = mod360( !radeg * cp_arr - !radeg * Transpose( temporary( rebin( cp_mean, n_bispect, nSim, /sample ) ) ) )
  cp_arr_meansub /= !radeg

  cp_cov = Transpose( cp_arr_meansub ) # cp_arr_meansub / nSim

  cp_err = fltarr( n_bispect )
  for i = 0, n_bispect - 1 DO cp_err[i] = Sqrt( cp_cov[i,i] )

  v2_mean = ABS( total( bispect_targ, 1 ) / total( bispect_cals, 1 ) )^2
  v2_arr_meansub = ABS( bs_arr )^2 - Transpose( temporary( rebin( v2_mean, n_bispect, nSim, /sample ) ) )
  
  v2_cov = Transpose( v2_arr_meansub ) # v2_arr_meansub / nSim
  v2_err = fltarr( n_bispect )
  for i = 0, n_bispect - 1 DO v2_err[i] = Sqrt( v2_cov[i,i] )

END



     
     
