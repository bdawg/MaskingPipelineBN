;; ----------------------------------------------
;;  This code uses a bootstrapping method to calculate the error
;;   and covariances between closure phases from bispectrum.
;;   Bootstrapping is a more accurate method at low signals to noise
;;   (errors tend to increase by ~ 50% using bootstrapping over a
;;   gaussian approximation).  
;;

  ;; The structure of our data taking is as follows:  We take sets of
  ;; images of a particularly object (target or calibrator) from which
  ;; we construct bispectrum.  These sets of data are kept array as
  ;; pointers in bispect_ptrs.

  ;; As part of our visual discrimination of data, we allow
  ;; each 'target set' to have associated with it particular
  ;; 'calibrator sets'; ie, not all calibrators are used on all
  ;; targets.

  ;; We wish to calculate the final 'science' bispectrum (
  ;; calibrated targets) using a bootstrapping method.  But this
  ;; requires combining various target sets of varying sizes, each of
  ;; which is calibrated with sets of data of various sizes.  We need
  ;; to take care in how we do this.

  ;; First, take look at the calibrators associated with each target
  ;; set, and we combine these by bootstrapping into a set of
  ;; nSim 'mean calibrators' for that target.  These are stored in
  ;; bootstrapped_cals, an array of size nSim x n_bispect x (target
  ;; set number).  

  ;; Among all the target sets, there are tot_targ bispectrum target
  ;; data, that each need to be calibrated (and which each have an
  ;; associated calibrator group to draw from).

  ;; The 'mean science bispectrum' would be the average of these
  ;; tot_targ bispectrum, each calibrated appropriately.
  
  ;; This is done by bootstrapping.  We do this many (n_large) times:
  ;; Draw tot_targ target bispectrum, calibrate each one with a randomly
  ;; selected 'mean calibrator' that was bootstrapped with its
  ;; respective calibrator.  The average of these tot_targ calibrated
  ;; target bispectrum is one statistically possible 'science bispectrum' from which
  ;; we build up the 'science bispectrum' distribution.  In this case,
  ;; n_large is arbitrarily chosen as n_targ_sets * nSim, so that each
  ;; bootstrap algorithm generates roughly the same number of bootstraps.

  ;; This sounds complex, but is only so because the numbers of targets and
  ;; calibrators can vary which requires a bit extra effort.  It
  ;; allows a very robust way to measure bispectrum means, errors, and
  ;; covariances without assumptions of gaussinity, which because
  ;; particularly important (very important) when signal to noise is
  ;; low or visibility amplitudes are low (for instance L-Dwarf Data).

;;  Author: David Bernat, dwb29@cornell.edu
;;

pro bootstrap_apm_stats, bispect_ptrs, cal4src, full_stats, bs_arr, NoCal=NoCal, nSim=nSim

  If( NOT keyword_set( nSim ) ) THEN nSim = 2000
  If( NOT keyword_set( NoCal ) ) THEN NoCal = 0

  n_sets = n_elements( bispect_ptrs )
  
  ;; ********************************************************
  ;; We need to know the number of closure phases in one data
  ;; measurement (usually 84)
  ;; *******************************************************
  first_bispect = *(bispect_ptrs[0])
  sz = size( first_bispect, /dimensions )
  n_bispect = sz[1]
  first_bispect = 0.

  ;; ***************************************************************
  ;;  This targ_2_grp maps from set index to target index.  It also
  ;;  will keep track of how many data points are taken for each set
  ;;  of target data, to properly weight their selection later.
  ;; ***************************************************************
  
  If( NoCal ) THEN BEGIN
     n_targ_sets = n_sets 
     For i = 0, n_sets - 1 DO BEGIN
        n_targ = (size( *(bispect_ptrs[i]), /dimensions ))[0]
        If( i EQ 0 ) THEN BEGIN
           targ_2_grp = replicate( i, n_targ )
        ENDIF ELSE BEGIN
           targ_2_grp = [ targ_2_grp, replicate( i, n_targ ) ] 
        ENDELSE
     ENDFOR
     ;; If there are calibrators
  ENDIF ELSE BEGIN
     n_targ_sets = 0
     For i = 0, n_sets - 1 DO BEGIN
        cals = where( cal4src[*,i] EQ 1, count )
        ;; If there's a calibrator with it, then its a target
        If( count EQ 0 ) THEN Continue

        n_targ = (size( *(bispect_ptrs[i]), /dimensions ))[0]
        If( n_targ_sets EQ 0 ) THEN BEGIN
           targ_2_grp = replicate( n_targ_sets, n_targ ) 
        ENDIF ELSE BEGIN
           targ_2_grp = [ targ_2_grp, replicate( n_targ_sets, n_targ ) ]
        ENDELSE
        n_targ_sets++
     ENDFOR
  ENDELSE 

  ;; ********************************************************
  ;;  From each target, take its associated calibrators, and bootstrap
  ;;  nSim 'mean calibrators' to calibrate each target with later.
  ;; ********************************************************

  Print, "Calculating Calibrator Statistics."
  bootstrapped_cals = complexarr( nSim, n_bispect, n_targ_sets )
  If( NoCal ) THEN BEGIN
     bootstrapped_cals = replicate( complex( 1.0, 0.0 ), NSim, n_bispect, n_targ_sets ) 
  ENDIF ELSE BEGIN
     i_targ_set = 0
     For i= 0, n_sets - 1 DO BEGIN
        cals = where( cal4src[*,i] EQ 1, count )
        ;; There are no calibrators for this data, don't analyse it.
        If( count EQ 0 ) THEN Continue
        ;; Otherwise, aggregate all its calibrators together
        For j = 0, n_elements( cals ) - 1 DO BEGIN
           If( j EQ 0 ) THEN BEGIN
              bispect_cals = *(bispect_ptrs[ cals[j] ]) 
           ENDIF ELSE BEGIN
              bispect_cals = Transpose([ [Transpose(bispect_cals)], [Transpose(*(bispect_ptrs[ cals[j] ]))] ])
           ENDELSE
        ENDFOR

        bootstrap_bispectrum, bispect_cals, boot_cals, nSim = nSim
        bootstrapped_cals[ *, *, i_targ_set ] = boot_cals
        i_targ_set++
     ENDFOR
  ENDELSE

  ;; ******************************************************
  ;;  From each target, create 'mean target' bispectrum
  ;;
  ;; ******************************************************

  Print, "Calculating Target Statistics."
  bootstrapped_targ = complexarr( nSim, n_bispect, n_targ_sets )
  
  If( NoCal ) THEN BEGIN
     For i = 0, n_sets - 1 DO BEGIN
        bispect_targ = *(bispect_ptrs[ i ])

        bootstrap_bispectrum, bispect_targ, boot_targ, nSim = nSim
        bootstrapped_targ[ *, *, i ] = boot_targ
     ENDFOR
  ENDIF ELSE BEGIN
     i_targ_set = 0
     For i= 0, n_sets - 1 DO BEGIN
        cals = where( cal4src[*,i] EQ 1, count )
        ;; There are no calibrators for this data, don't analyse it.
        If( count EQ 0 ) THEN Continue
        ;; Otherwise, it is a target.  Bootstrap it.
        bispect_targ = *(bispect_ptrs[ i ]) 

        bootstrap_bispectrum, bispect_targ, boot_targ, nSim = nSim
        bootstrapped_targ[ *, *, i_targ_set ] = boot_targ
        i_targ_set++
     ENDFOR
  ENDELSE

  ;; ******************************************************
  ;;  Now bootstrap the 'mean science bootstraps'.  To generate one
  ;;  mean science bispectrum, draw one 'mean target bispectrum',
  ;;  calibrate it with one of its generated 'mean
  ;;  calibrators'.  That's a mean science bispectrum.  We want to
  ;;  weight each of the target sets by the number of data points in
  ;;  each set.  (An alternative is to weigh them by the inverse of
  ;;  their standard deviations, but this is more straightforward.)
  ;;
  ;;  This code bootstrapped n_targ_sets * nSim total science bispectrum
  ;; *********************************************************

  ;; The code before is similiar to that in bootstrap_bispect.pro, but
  ;; bootstrap_bispect is a general bootstrapper.
  n_final_sim = nSim * n_targ_sets
  tot_targ = n_elements( targ_2_grp )

  bs_arr = complexarr( n_final_sim, n_bispect )
  print, "Starting Bootstrapping Simulation of Calibrated Targets..."

  ;; Because of a pecularity of how IDL handles array subsets, the
  ;; bootstrapped targs and cals (3D array) need to be reformed to 2D
  ;; arrays.
  boot_cals_rf = complexarr( nSim*n_targ_sets, n_bispect )
  For i = 0, n_targ_sets - 1 DO boot_cals_rf[ (i*nSim):((i+1)*nSim-1), * ] = bootstrapped_cals[*,*,i]
  bootstrapped_cals = 0.
  boot_targ_rf = complexarr( nSim*n_targ_sets, n_bispect )
  For i = 0, n_targ_sets - 1 DO boot_targ_rf[ (i*nSim):((i+1)*nSim-1), * ] = bootstrapped_targ[*,*,i]
  bootstrapped_targ = 0.

  rn_set  = floor( randomu( seed, n_final_sim ) * tot_targ )
  rn_targ = floor( randomu( seed, n_final_sim ) * nSim )
  rn_cals = floor( randomu( seed, n_final_sim ) * nSim ) 
  rn_cals = targ_2_grp[ rn_set ] * nSim + rn_cals
  rn_targ = targ_2_grp[ rn_set ] * nSim + rn_targ

  bs_arr = boot_targ_rf[ rn_targ, * ] / boot_cals_rf[ rn_cals, * ]


  ;; *****************************
  ;;  Statistics
  ;; *****************************

  ;; Calculate the statistics of:
  ;;  1) Science bispectrum
  ;;  2) Aggregated Calibrators for each Target
  ;;  3) Uncalibrated Targets
  ;;  4) Each Calibrated Target

  ;; NOTE: Swaping #2 for Stats of Each Calibrator set, ie,
  ;; non-aggregated, would probably be more illuminating.  No time
  ;; right now.

  ;; First, rearrange the bootstraps back to a readable style
  bootstrapped_cals = complexarr( nSim, n_bispect, n_targ_sets )
  For i = 0, n_targ_sets - 1 DO bootstrapped_cals[ *, *, i ] = boot_cals_rf[ (i*nSim):((i+1)*nSim-1),*]
  boot_cals_rf = 0.
  bootstrapped_targ = complexarr( nSim, n_bispect, n_targ_sets )
  For i = 0, n_targ_sets - 1 DO bootstrapped_targ[ *, *, i ] = boot_targ_rf[ (i*nSim):((i+1)*nSim-1),*]
  boot_targ_rf = 0.

  ;; Second, get a 'blank' statistics structure
  
  blank_bs = complexarr( 2, 84 ) + 1.0 
  crunch_cp_stats, blank_bs, blank_stats

  ;; Calibrator Stats
  cals_stats = replicate( blank_stats, n_targ_sets )
  
  For i = 0, n_targ_sets - 1 DO BEGIN
     crunch_cp_stats, bootstrapped_cals[*,*,i], cs
     cals_stats[i] = cs
  ENDFOR

  ;; Uncalibrated Target Stats
  uncal_targ_stats = replicate( blank_stats, n_targ_sets )
 
  For i = 0, n_targ_sets - 1 DO BEGIN
     crunch_cp_stats, bootstrapped_targ[*,*,i], cs
     uncal_targ_stats[i] = cs
  ENDFOR

  ;; Each Target Set, Calibrated
  ;; The 'target' bispectrum still need to be calibrated.  Since both
  ;; the targ and cal each have nSim bootstraps, and the bootstraps
  ;; are independent (boot 0 of targ is not correlated with boot 0 of
  ;; cal), we can take the easy way out and divide the two sets
  ;; without sacrificing any randomness.
  targ_stats = replicate( blank_stats, n_targ_sets )

  For i = 0, n_targ_sets - 1 DO BEGIN
     crunch_cp_stats, bootstrapped_targ[*,*,i] / bootstrapped_cals[*,*,i], cs
     targ_stats[i] = cs
  ENDFOR

  ;; Science Stats
  crunch_cp_stats, bs_arr, science_stats

  full_stats = { cals_stats: cals_stats, targ_stats: targ_stats, uncal_targ_stats: uncal_targ_stats, science_stats:science_stats } 

END



  
  
