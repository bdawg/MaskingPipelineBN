function generate_noise_bs, cp_mean, cp_cov, nSim=nSim, seed=seed
  
  If( NOT keyword_set( nSim ) ) THEN nSim = 10000
  If( NOT keyword_set( seed ) ) THEN seed = 2L
 
  n_bispect = n_elements( cp_mean )

  ;; Generate Noise matching Covariances
  Cov_EI = EigenQL( cp_cov, Eigenvectors = Cov_evecs, /double )
  sig_z  = randomn( seed, nSim, n_bispect ) # diag_matrix( sqrt( cov_ei ) )
  sig_cp = sig_z # Transpose( Cov_evecs )

  ;;sig_cp += cp_mean
  return, complex( Cos( sig_cp / !radeg ), Sin( sig_cp / !radeg ) )
END
