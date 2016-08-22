pro bootstrap_bispectrum, bispect, boot_arr, nSim=nSim, seed=seed
  
  If( NOT keyword_set( nSim ) ) THEN nSim = 2000
  ;; If seed is not set, leave it undefined.

  sz = size( bispect, /dimensions )
  n_data    = sz[0]             ; Number of data points
  n_bispect = sz[1]          
  ;; This is a semi-arbitrary number that I pulled from my own anecdotal
  ;; testing.  This is the number of elements in the operational array.
  ;; Too large and it crushes the RAM. (1 float = 4 bytes)
  limit = 45.e6

  n_sim_per = Ceil( limit / n_data  / n_bispect ) < nSim
  n_bunches = Ceil( 1.0*nSim / n_sim_per )

  boot_arr = complexarr( nSim, n_bispect )
  print, "Starting Bootstrapping Simulation..."
  done_boots = 0
  For n = 0, n_bunches - 1 DO BEGIN
     ;; We don't want to generate more than nSim bootstraps.
     n_sim_this = ( nSim - done_boots ) < n_sim_per

     ;; Seed is updated with each call, so this is does generate new
     ;; pseudo-random numbers each time.  Also is seed is undefined.
     rn = floor( randomu( seed, n_sim_this * n_data ) * n_data )

     ;; This used to be done on a per-bispectrum level -- otherwise
     ;; the array sizes may be large enough to crash the RAM.  But
     ;;  now the bispectrum arrays are on the order of a few hundred, not a
     ;;  few ten thousand.

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

     boot_arr[(done_boots):(done_boots+n_sim_this-1),*] = total( reform( temporary( bispect[rn,*] ), n_sim_this, n_data, n_bispect, /overwrite ), 2 ) / n_data
     done_boots += n_sim_this
     print, "Finished Bunch ", n+1, " of ", n_bunches, " bunches.  "
     print, "Finished ", done_boots, " bootstraps." 
  ENDFOR
END






