;; --------------------------------------------------
;;
;;  This script calculates Detection Contrast Limits for t3data.
;;
;;  This code approximates the following:  Take your closure phases,
;;fit for the best binary, and compare it to the F distribution
;;calculated from emp_f_dist.pro.
;;
;;  With what confidence is the binary Phi( rho, theta, r )
;;detectable?  Well, for each nSim set of Noise_CPs which represent
;;single star measurements, we generate noisy binary measurements:
;;
;; Bin_CP = Noise_CP + Phi( rho, theta, r )
;;
;; Then, 'fit' that binary with the known parameters, calculate the F
;;of the fit, convert that to significance, and average that
;;significance over all nSim set of noisy CPs.  That's the
;;confidence of detecting those parameters.  
;;  The 99.5% detection limit is then the contrast ratio dimmest
;;contrast with an average confidence of .995.
;;
;; Keywords:
;;
;;  NOISE_BS: Same noise_bs as passed into emp_f_dist.pro
;; 
;;  CP_COV or DATALIKES: One or the other, or neither.  Datalikes are
;;expected to be mean zero.  That is, the maximum of the likelihood
;;should be shifted to zero.
;;
;;  BEST_GRID_F: Retrieved from emp_f_dist.pro output.
;;  NOISE_NULLCHI2: Retrieved from emp_f_dist.pro. (See comments below.)
;;
;; 
;; Usage Notes:
;;
;;  NSIM is SQUARE: See the Speed Trick comment below.  This routine
;;runs significantly faster is nSim is a square number.  ie., nSim = 100^2
;;
;;  WHY BEST_GRID_F?: This routine doesn't do a gradient search
;;to find best parameters, so the finished F distribution from
;;emp_f_dist.pro isn't the proper distribution to use.  This
;;routine sticks to parameters on the parm grid.  To keep things
;;consistent, the empifical F distribution is built from the best F
;;values on the parameter grid.  Inspection of emp_f_dist shows that
;;the output Best_Grid_F is exactly this distrbution.
;;
;;
;; Desirable Improvements:
;;
;;  Dynamic Crat: It'd be nice to automatically test with
;;contrast ratios down to the 99.5% limit, rather than be locked into
;;the range passed in by ParmGrid.
;;
;;  Explicit Parm Grid: Rather than passing in a min-max and number,
;;pass in the explicit parameters values to calculate on.  ie., pass
;;in sep = [ 45.0, 65.0, 125.0, 200.0, 300.0, 400.0 ]
;;
;; ---------------------------------------------------------------------

pro emp_crat_limits, t3data, parmlimits, Noise_bs, cp_cov=cp_cov, datalikes=datalikes, Best_Grid_F, Noise_nullchi2, ret_struct

  sz = size( Noise_bs, /dimension )
  nSim = sz[0]
  n_bispect = sz[1]

  ;; Template of UV-data passed into T3Models
  t3template = t3data

  ;; Prepare Data.  Convert from BS to CP.  Noise_Z is statistically independent and
  ;; what is fit.  We do this so we don't need to invert the
  ;; covariance matrix.
  noise_cp = atan( noise_bs, /phase ) * !radeg
  noise_z  = noise_cp  ;; Same thing if not using cov
  Z_err    = t3data.t3phierr

  ;; -------------------------
  ;; Grid Parameters
  ;; -------------------------

  minSep = ParmLimits[0].min
  maxSep = ParmLimits[0].max
  nSep   = ParmLimits[0].n

  minAz = 0   ;;ParmLimits[1].min
  maxAz = 360 ;;ParmLimits[1].max
  nAz   = 12  ;;ParmLimits[1].n

  minCRat = ParmLimits[2].min
  maxCRat = ParmLimits[2].max
  nCRat   = ParmLimits[2].n

  Sep  = minSep + (maxSep - minSep) * findgen(nSep)/nSep
  Az   = minAz  + (maxAz - minAz) * findgen(nAz) / nAz
  Crat = 10.^( minCrat + (maxCrat - minCrat) * findgen(nCrat)/nCrat )

  Sep = [ 65., 85., 105., 125., 145., 165., 185., 225., 265., 305., 345., 385., 425. ]
  nSep = n_elements( Sep )
  minSep = min( Sep )
  maxSep = max( Sep )

  ;; Don't bother doing the equal contrast cases
  If( minCrat EQ 0.0 ) THEN BEGIN
     minCrat = Crat[ 1 ]
     Crat = Crat[1:(nCrat-1)]
     nCrat = nCrat - 1
  ENDIF

  ;; --------------------------------------------
  ;; Using Covariance matrix or Data likelihoods?
  ;; --------------------------------------------

  use_angle_diff = 1.0  
  If( Keyword_Set( CP_COV ) ) THEN BEGIN
     use_angle_diff = 0.0
     Cov_EI = EigenQL( cp_Cov, Eigenvectors = Cov_evecs, /double )
     Noise_Z = Noise_CP # Cov_Evecs
     CP2Z_MTX = Transpose( Cov_Evecs ) ;; Matrix to go from CP to Z for models.
     Z_err = Sqrt( Cov_EI )
  ENDIF

  IF( Keyword_set( Datalikes ) ) THEN BEGIN
     use_angle_diff = 0.0
  ENDIF

  ;; Create variables expected to exist
  If( NOT Keyword_Set( CP_Cov ) ) THEN BEGIN
     cp_cov = 0.0
     CP2Z_MTX = 0.0
  ENDIF

  If( NOT Keyword_set( datalikes ) ) THEN datalikes = 0.0

  ;; ---------------
  ;;  Speed Trick
  ;; ---------------

  ;; It is faster to search through a sorted 100x100 array than a
  ;; 10000x1 array.  Also, referencing an array by X[0,*] takes a 
  ;; significant amount of time.

  F_dist = Best_grid_F[ sort( Best_Grid_F ) ]
  If( Sqrt( nSim ) EQ Floor( Sqrt( nSim ) ) ) THEN BEGIN
     n_side = Sqrt( nSim )
     F_dist_square = reform( F_dist, n_side, n_side )
     F_dist_0 = F_dist_square[0,*]
     square = 1.0
  ENDIF ELSE square = 0.0


  ;; Output Variable
  psig_cube = fltarr( nSep, nAz, nCrat )

  tb = systime(1)
  For i = 0, nSep - 1 DO BEGIN
     For j = 0, nCrat - 1 DO BEGIN
        For k = 0, nAz - 1 DO BEGIN

           ;; For each Simulated CP, calculate the F statistic for
           ;; this model, and the probability of significance.
           p0 = [ Sep[i], Az[k], Crat[j] ]
           
           t3model = binary_t3data_fast( p0, t3data=t3template )

           ;; Here's what we do: (except for Datalikes)
           ;; Sim_CP is a mock observation of only noise.
           ;; Bin_CP = Model_CP + Sim_CP is a mock binary.
           ;; Then, we calculate F = (Null - Best) / Best

           ;; We assume the Best Fit is the True Fit.
           ;;  In that case, note:
           ;; Best Chi2 = Sum[ (Model_CP + Sim_CP - Model_Cp)^2 ] 
           ;;           = Sum[ Sim_CP^2 ] = What was called Null Chi2
           ;;           in emp_f_dist.  Instead of Calculating this again,
           ;;           we passed it in to this function.
           
           ;; Null Chi2 = Sum[ Bin_CP^2 ] = Sum[ (Model_CP + Sim_CP)^2 ]
           
           ;; Finally, note:
           ;; Null - Best = 
           ;;  Sum[ Model_CP^2 ] + 2 * Sum [ Model_CP * Sim_CP ]

           ;; We calculate this quantity instead for the numerator.

           If( Keyword_Set( datalikes ) ) THEN BEGIN
              ;; DataLikes
              res2_arr = fltarr( n_bispect, nSim )
              For q = 0, n_bispect - 1 DO res2_arr[q,*] = -alog( interpol( (*datalikes[q]).l_vec, (*datalikes[q]).m_vec, -mod360( Noise_Z[*,q] + (t3model.t3phi)[q] ) ) > 1.D-300 ) ;; This is correct.
              chi2_arr = total( res2_arr, 1 )
              dChi2 = chi2_arr - Noise_NullChi2              
           ENDIF ELSE IF( KeyWord_Set( CP_COV ) ) THEN BEGIN
              ;; Covariances
              Model_Z =  CP2Z_MTX # (t3model.t3phi)
              
              Model_bit = Model_Z / Z_err^2 ;; Useful
              Term2 = replicate( Total( (Model_Z/Z_err)^2 ), nSim )

              Model_bit = Transpose( Rebin( Model_bit, n_bispect, nSim ) ) ;; Arrayify it.
              Term1 = 2 * Total( Model_bit * Noise_Z, 2 )
              dChi2 = Term1 + Term2
           ENDIF ELSE BEGIN
              ;; Errors
              Model_bit = (t3model.t3phi) / Z_err^2 ;; Useful
              Term2 = replicate( Total( ((t3model.t3phi)/Z_err)^2 ), nSim )

              Model_bit = Transpose( Rebin( Model_bit, n_bispect, nSim ) ) ;; Arrayify it.
              Term1 = 2 * Total( Model_bit * Noise_Z, 2 )
              dChi2 = Term1 + Term2
           ENDELSE

           F = dChi2 / Noise_nullChi2

           ;; --------------------------------------
           ;; Convert these F Stats to Signifcances
           ;; -------------------------------------
           psig = fltarr( nSim )
           If( square NE 0.0 ) THEN BEGIN
              For l = 0, nSim - 1 DO BEGIN
                 i1 = 0 > ( total( F[l] GT F_dist_0 ) - 1 )
                 i2 = total( F[l] GT F_dist_square[*,i1] )
                 psig[l] = n_side * i1 + i2 
              ENDFOR
           ENDIF ELSE For l = 0, nSim - 1 DO psig[l] = total( F[l] GT F_dist )

           ;; -------------------------------------
           ;;  What we want is the avg psig over all nSim
           ;; -------------------------------------
           psig_avg = 1.0 * total( psig ) / nSim^2
           psig_cube[ i, k, j ] = psig_avg
        ENDFOR
     ENDFOR
     te = systime(1)
     print, "Simulation running:", 100.*(i+1)/nSep, "% complete in ", te-tb, " seconds."
     print, "Simulation will be complete at: ", systime( 0, tb + (te-tb)*nSep/(i+1) )
  ENDFOR

  ;;   ***** This section is now calculated on an as needed basis elsewhere in the code ****
  ;;   Now, we want an achievable contrast limit for a given separation.
  ;;   So, we average over the thetas.
  ;;   Then, we extrapolate to what contrast limit will give us our desired
  ;;   confidence limit.

  ;;psig_avg = total( psig_cube, 2 ) / nAz

  ;;Crat_Limits = fltarr( nSep )
  ;;Mag_Limits  = fltarr( nSep )
  ;;For i = 0, nSep - 1 DO BEGIN
  ;;   Crat_lev = interpol( Crat, psig_avg[i,*], conf_limit )

  ;; ind_Crat is now the interpolated index of the Crat vector
  ;;   Crat_Limits[i] = Crat_lev
  ;;   Mag_limits[i]  = 2.5 * alog( Crat_lev ) / alog( 10. )
  ;; ENDFOR
  ;;  *************** END *******************

  ;; ---------------------------------------------
  ;;  Finally, convert the parm grid into a format to store along with
  ;;  the results.
  ;; --------------------------------------------
  parm_cube = fltarr( 3,nSep, nAz, nCrat)
  For i = 0, nSep - 1 DO For j = 0, nAz - 1 DO FOR k = 0, nCrat - 1 DO $
     parm_cube[*,i, j, k] = [ Sep[i], Az[j], Crat[k] ]

  ret_struct = { parm_cube: parm_cube,  psig_cube: psig_cube }
END
