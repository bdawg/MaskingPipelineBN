;; --------------------------------------------------------
;;
;;  The purpose of this code is to generate the false alarm
;;probability of a binary detection.
;;
;;  This code finds the best fitting models to each of the sets
;;of closure phases/bispectrum passed in through the noise_bs array.
;;
;;  NOISE_BS: Mock observations of a single star with the same noise 
;; properties as t3data.  Each mock observation is fit and delta-Chi2
;; and % improvement are recorded (along with the best fit
;; parameters).  This generates a distribution of these metrics that
;; arise from single star observations.  By comparing the delta-Chi2 or
;; % improvement to this distribution, one gets the probability of
;; false alarm; i.e., that the delta-Chi2 value arose not from a binary
;; but from a single star.
;;
;;  This code works with or without a Covariance matrix (CP_COV) and
;;with data likelihoods instead of the t3data (if DATALIKES are passed in).
;;
;;  This is called 'EMP_F_DIST' because % improvement of Chi2 is an F statistics.
;;
;;  Keywords:
;;
;;   DYNAMIC_CRAT_OFF: The best fitting binary model has a contrast
;;ratio r that is approximately equal to the closure phases values in
;;radians.  With this keyword set, the grid of contrast ratios
;;searches is dynamically generated based on the set of closure phases
;;being fit.  This ensures that the grid search searches the most
;;likely best fits.  Setting this keyword to 0 is RECOMMENDED. (i.e,
;;ON).
;;
;;    DISPLAY_GRAPH: Display a graph of every 10th Best-Fit, Chi2
;;Contour, and current Histogram of % Improvement Values.  Useful for
;;making sure everything is working properly.
;;
;;  Usage Notes:
;;
;;   DATALIKES: These are expected to have their mean subtracted off.
;;In other words, the max of the likelihood curve should be at ZERO.
;;
;;
;;  Author: David Bernat, dbernat@physics.cornell.edu
;;
;; -------------------------------------------------------------   

pro emp_f_dist, t3data, parmlimits, noise_bs, ret_struct, CP_Cov=CP_Cov, datalikes=datalikes, Dynamic_Crat_Off=Dynamic_Crat_Off, DisplayGraph=DisplayGraph

  If( NOT Keyword_Set( DisplayGraph ) )     THEN DisplayGraph = 0
  If( NOT keyword_set( Dynamic_Crat_Off ) ) THEN Dynamic_Crat_Off = 0

  If( Keyword_Set( CP_Cov ) AND KeyWord_Set( Datalikes ) ) THEN BEGIN
     message, "One or the other, pal.  CP COvariances or Data Likelihoods, but not both."
  ENDIF

  sz = Size( noise_bs, /dimension )
  nSim = sz[0]
  n_bispect = sz[1]
  
  ;; Output Variable
  BestFitOutput = fltarr( nSim, 7 )

  ;; Template of UV-data passed into T3Models
  t3template = t3data

  ;; Prepare Data.  Convert from BS to CP.  Noise_Z is statistically independent and
  ;; what is fit.  We do this so we don't need to invert the
  ;; covariance matrix.
  noise_cp = atan( noise_bs, /phase ) * !radeg
  noise_z  = noise_cp  ;; Same thing if not using cov
  t3sim  = { t3phi: fltarr( n_bispect ), t3phierr: (t3data.t3phierr ) }
  t3null = t3sim

  ;; Using Covariance matrix or Data likelihoods?
  use_angle_diff = 1.0  
  If( Keyword_Set( CP_COV ) ) THEN BEGIN
     use_angle_diff = 0.0
     Cov_EI = EigenQL( cp_Cov, Eigenvectors = Cov_evecs, /double )
     Noise_Z = Noise_CP # Cov_Evecs
     CP2Z_MTX = Transpose( Cov_Evecs ) ;; Matrix to go from CP to Z for models.
     t3sim.t3phierr = Sqrt( Cov_EI )
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

;; -----------------------------------------------------------------------------
;;                           F-Fit for A BINARY SYSTEM
;; -----------------------------------------------------------------------------

  minSep = ParmLimits[0].min
  maxSep = ParmLimits[0].max
  nSep   = ParmLimits[0].n

  minAz = ParmLimits[1].min
  maxAz = ParmLimits[1].max
  nAz   = ParmLimits[1].n

  nCRat   = ParmLimits[2].n

  Sep  = minSep + (maxSep - minSep) * findgen(nSep)/nSep
  Az   = minAz  + (maxAz - minAz) * findgen(nAz) / nAz

  If( Dynamic_Crat_Off ) THEN BEGIN
     minCRat = ParmLimits[2].min
     maxCRat = ParmLimits[2].max

     Crat = 10.^( minCrat + (maxCrat - minCrat) * findgen(nCrat)/nCrat )
  ENDIF

;; ---------------------------------------
;;   Prep for the grid searches 
;; --------------------------------------

;; ---- constraints on the parameters in MPFIT ----

  pi = replicate({fixed:0, limited:[0,0], limits:[0.,0.], relstep:.01}, 3)
  
  ;; Limit PA to within Range of Opti's

  pi[1].limited[0] = 1
  pi[1].limited[1] = 1
  pi[1].limits[0]  = Az[ 0 ]
  pi[1].limits[1]  = Az[ (nAz-1) ]

  pi[0].limited[0] = 1
  pi[0].limited[1] = 1
  pi[0].limits[0]  = Sep[ 0 ]
  pi[0].limits[1]  = Sep[ (nSep-1) ]

  If( Dynamic_Crat_Off ) THEN BEGIN
     pi[2].limited[0] = 1
     pi[2].limited[1] = 1
     pi[2].limits[0]  = Crat[ 0 ]
     pi[2].limits[1]  = Crat[ (nCrat-1) ]
  ENDIF

  tbegin = systime(1)
  For n = 0, nSim - 1 DO BEGIN
     t3sim.t3phi = Reform( Noise_Z[n,*] ) ;; Used for Chi2 Calc

     ;; The Datalikes have to be shifted so their max is at the mock
     ;; value.  They are shifted back at the end.
     IF( Keyword_Set( datalikes ) ) THEN BEGIN
        For i = 0, n_bispect - 1 DO BEGIN
           m_vec = (*datalikes[i]).m_vec
           l_vec = (*datalikes[i]).l_vec
           m_vec = mod360( m_vec + (t3sim.t3phi)[i] )
           s = sort( m_vec )
           m_vec = m_vec[s]
           l_vec = l_vec[s]
           *datalikes[i] = { m_vec: m_vec, l_vec: l_vec }
        ENDFOR
     ENDIF
     

     ;; ****************************************
     ;;  Calculate Fit Range for Crat Dynamically
     ;; ****************************************

     If( NOT Dynamic_Crat_Off ) THEN BEGIN
        CP_sim = Reform( Noise_CP[n,*] ) 
        
        ;; For high Crat, CP ~ 1/Crat.
        ;; For Crat ~ (1+ep), CP ~ 1/ep.
        ;; From these the below follows.

        ;; So look for the 10th and 74th (of 84) bounded simulated CPs and use
        ;; them as an estimate for 1/Crat.  Piece together a high Crat
        ;; and a Crat ~ 1 search space.

        t3abs = ABS( cp_sim ) 
        s = sort( t3abs )
        
        ;; Low SNR can have huge CPs, driving hi_crat_min below 1.0.
        ;; Firstly, Crat is always greater than 1.  Secondly, it would
        ;; mean that the 'high crat' range is dipping into/overlapping
        ;; with the 'low crat' range.

        lowind = ceil( .11 * n_bispect ) ;; ~ 10 for 84 bispect
        hi_crat_min = ( !radeg / t3abs[s[n_bispect-lowind]] )
        hi_crat_max = !radeg / t3abs[s[lowind]]

        low_crat_min = hi_crat_max / (hi_crat_max - 1. ) 
        low_Crat_max = hi_crat_min / (hi_crat_min - 1. )
        
        If( ( hi_crat_min LT 1.0 ) OR ( hi_crat_min LT low_crat_max ) ) THEN BEGIN
           ;; Low SNR Region
           minCrat = 1e-4                       ;; log scale
           maxCrat = alog10( hi_crat_max ) + .4 ;; A fudge factor 
           
           Crat = 10.^( minCrat + (maxCrat - minCrat) * findgen(nCrat)/nCrat )
        ENDIF ELSE BEGIN
           ;; High SNR Region
           hi_crat_min = alog10( hi_crat_min )
           hi_crat_max = alog10( hi_crat_max )
           low_crat_min = alog10( low_crat_min )
           low_crat_max = alog10( low_crat_max )
           dum1 = Floor( nCrat / 5. )
           dum2 = nCrat - dum1
           C1 = 10.^( low_crat_min + (low_crat_max - low_crat_min) * findgen(dum1)/dum1 )
           C2 = 10.^( hi_crat_min + (hi_crat_max - hi_crat_min) * findgen(dum2)/dum2 )
           
           ;; This is the total Crat space searched
           Crat = [ C1, C2 ] 
           minCrat = min( Crat )
           maxCrat = max( Crat )
        ENDELSE

        ;; Set the correct MPFIT limits
        pi[2].limited[0] = 1
        pi[2].limited[1] = 1
        pi[2].limits[0]  = Crat[ 0 ]
        pi[2].limits[1]  = Crat[ (nCrat-1) ]
     ENDIF

     ;; If a given location has chi2 > nullchi2, don't bother searching.
     nullchi2 = binary_chi2_fast( [ 1.0, 1.0, 1.0 ], t3model=t3null, t3data=t3sim, datalikes=datalikes, use_angle_diff=use_angle_diff )
     
     ;; This is some array-fu that calculates the chi2 slightly faster
     ;; by taking advantage of the IDL speed of operations on arrays
     n_d = nSep*nAz*nCrat
     parmgrid = fltarr( 3, n_d )
     For i = 0, nSep - 1 DO For j = 0, nAz - 1 DO For k = 0, nCrat - 1 DO $
        parmgrid[ *, i + nSep*j + nSep*nAz*k ] = [ Sep[i], Az[j], Crat[k] ] 
     
     t3modelCPs_arr = binary_t3data_many( parmgrid, t3data=t3template, /Return_CPs_Only )
     
     If( Keyword_Set( CP_COV ) ) THEN BEGIN
        ;; Covariances
        t3modelZs_arr = CP2Z_MTX # ( t3modelCPs_arr )
        t3simZs_arr   = Rebin( t3sim.t3phi, n_bispect, n_d, /sample )
        t3simZerr_arr = Rebin( t3sim.t3phierr, n_bispect, n_d, /sample )
        chi2_arr = Total( ( (t3modelZs_arr - t3simZs_arr) / t3simZerr_arr )^2, 1 )
     ENDIF ELSE If( Keyword_Set( Datalikes ) ) THEN BEGIN
        ;; Data Likelihoods
        res2_arr = fltarr( n_bispect, n_d )
        For i = 0, n_bispect - 1 DO res2_arr[i,*] = -alog( interpol( (*datalikes[i]).l_vec, (*datalikes[i]).m_vec, t3modelCPs_arr[i,*] ) > 1D-300 )
        chi2_arr = total( res2_arr, 1 )
     ENDIF ELSE BEGIN
        ;; Straight Chi2 
        t3sim_arr    = Rebin( t3sim.t3phi, n_bispect, n_d, /sample )
        t3simerr_arr = Rebin( t3sim.t3phierr, n_bispect, n_d, /sample )
        chi2_arr = Total( ( (t3modelCPs_arr - t3sim_arr) / t3simerr_arr )^2, 1 )
     ENDELSE

     ;; Search only the lowest Chi2 in the bunch
     firstchi2 = Min( Chi2_arr, ind )
     p0 = parmgrid[ *, Ind ]
     
     ;; ---------------------
     ;;  Begin Gradient Search
     ;; ---------------------

     ;; --- Additional Arguments for mp_binary_func_fast ---
     fa = {t3data:t3sim, evecs: CP2Z_MTX, t3template:t3template, datalikes: datalikes, use_angle_diff: use_angle_diff }
     
     ;; Gradient Search
     pout = mpfit( 'mp_binary_func_fast', p0, functargs=fa, parinfo=pi, status=status, errmsg=errmsg, quiet = 1 )
     
     If( status LE 0 ) then message, errmsg
     
     ;; Get real chi2: mpfit has a bug that returns the correct
     ;; params, but wrong chi2. (Wrong Chi2 from next to last best fit.)
     chi2 = total( ( call_function( 'mp_binary_func_fast', pout, _EXTRA=fa ) )^2.0 )
     
     ;; See Bevington-Robinson for information on F Statistics.
     ;;  In practice, we scale the errors so that the Best Fit Reduced
     ;;  Chi2 = 1.  In that case, the calculated 
     ;;  DeltaChi2 = (NullChi2 - Chi2 ) / Chi2 * DOF = F * DOF
     ;;  Since DOF is constant, we can compare directly the F of the fitted
     ;;  data to the F distribution generated here.

     F = ( nullchi2 - chi2 ) / chi2
     BestFitOutput[n,*] = [ pout[0], pout[1], pout[2], nullchi2, firstchi2, chi2, F ]
     tend = tbegin + (systime(1) - tbegin) * NSim / ( n + 1 )
     Print, "Simulation ", (n+1), " is finished.  Expect to finish all ", NSim, " simulations at: ", systime( 0, tend ), ".  Avg Time: ", (systime(1) - tbegin)/(n+1) 
     print, "Simulation ", n, ": ", Transpose( BestFitOutput[n,*] )

     ;; Display Nifty Graphics
     If( Keyword_Set( Display_Graph ) AND (n mod 10 EQ 0 ) ) THEN BEGIN
        !p.multi = [ 0, 1, 3 ]

        ;; Current Histogram
        bin_per_int = 100.
        x = BestFitOutput[0:n,6] * bin_per_int
        h = histogram( x )
        b = findgen( n_elements( h ) ) + min( x ) 
        plot, b / bin_per_int, h / (n+1)

        ;; Current Best Fit
        If( KeyWord_Set( Datalikes ) ) THEN BEGIN
           ;; Datalikes
           ;; Contour Plot of each likelihood curve
           log_cont = fltarr( n_bispect, n_elements( (*datalikes[0]).m_vec ) )
           For i = 0, n_bispect - 1 DO log_cont[i,*] = sqrt( -alog( (*datalikes[i]).l_vec > 1d-300) ) < 5.   
           contour, log_cont, findgen(n_bispect), (*datalikes[0]).m_vec, /fill, nlevel=5, yrange=[-40,40]
           t3best = binary_t3data_fast( pout, t3data=t3template )
           oplot, (t3best.t3phi)[s], color=65000
           oplot, (t3best.t3phi)[s], color=65000, /psym
           oplot, replicate( 0.0, 84 ), color=250
        ENDIF ELSE IF ( Keyword_set( CP_COV ) ) THEN BEGIN
           ;; Covariances
           t3best = binary_t3data_fast( pout, t3data=t3template )
           ploterr, t3sim.t3phi, t3sim.t3phierr, color=250
           oplot, CP2Z_MTX # t3best.t3phi, /psym, color=65000
           oplot, fltarr( n_bispect ), color=250
        ENDIF ELSE BEGIN
           ;; T3Data Errors
           t3best = binary_t3data_fast( pout, t3data=t3template )
           ploterr, Noise_CP[n,*] , t3sim.t3phierr, color=250
           oplot, t3best.t3phi, /psym, color=65000
           oplot, fltarr( n_bispect ), color=250
        ENDELSE

        ;; Chi2 Contours
        minall = fltarr(nSep,nAz)
        for i = 0, nSep-1 do begin
           for j = 0, nAz - 1 do begin
              minall[i,j] = MIN( (reform( chi2_arr, nSep, nAz, nCrat ))[i,j,*] ) < nullchi2
              minall[i,j] = ( minall[i,j] - nullchi2 ) / nullchi2
           ENDFOR
        ENDFOR
        
        contour, minall, Sep, Az, /fill, tit=title , nlevel=20, xtitle='angular sep (mas)', ytitle = 'azimut (degrees)', xrange=[minSep, maxSep], yrange=[minAz, maxAz]
        contour, minall, Sep, Az, tit='F', nlevel=5, /overplot, xrange=[minSep, maxSep], yrange=[minAz, maxAz]
        plots, pout[0], pout[1], psym=2, color=65000, symsize=3
        
     ENDIF
     ;; END Nifty Display
     
     ;; Shift Data Likelihoods Back to Mean Zero
     IF( Keyword_Set( datalikes ) ) THEN BEGIN
        For i = 0, n_bispect - 1 DO BEGIN
           m_vec = (*datalikes[i]).m_vec
           l_vec = (*datalikes[i]).l_vec
           m_vec = mod360( m_vec - (t3sim.t3phi)[i] )
           s = sort( m_vec )
           m_vec = m_vec[s]
           l_vec = l_vec[s]
           *datalikes[i] = { m_vec: m_vec, l_vec: l_vec }
        ENDFOR
     ENDIF

     ;; Move on to next set of mock noise.
  ENDFOR

  Print, "DONE!!!  Phew."
  ret_struct = { BestFitOutput: BestFitOutput, noise_bs: noise_bs, cp_cov: cp_cov, datalikes: datalikes }

end
