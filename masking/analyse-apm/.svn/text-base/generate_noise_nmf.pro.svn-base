function generate_noise_nmf, datalikes, means=means, sys_t3=sys_t3, nSim=nSim, seed=seed
  
  If( NOT keyword_set( nSim ) )  THEN nSim = 10000
  If( NOT keyword_set( seed ) )  THEN seed = 2L
  
  ;; Assumes datalikes are mean-zero already
  If( NOT keyword_set( means ) ) THEN means = fltarr( n_elements( datalikes ) )

  ;; Noise CP = signal noise + systematic signal + systematic noise
  
  n_bispect = n_elements( datalikes )
  ;; Generate Signal Noise (i.e., measurement noise)
  rnds = randomu( seed, nSim, n_bispect )
  mock_cp = rnds * 0.0
  For i = 0, n_bispect - 1 DO BEGIN
     m_vec = mod360( (*datalikes[i]).m_vec - means[i] )
     s = sort( m_vec )
     m_vec = m_vec[s]
     l_vec = (*datalikes[i]).l_vec[s]

     c_vec = l_vec * 0.0
     For j = 1, n_elements( c_vec ) - 1 DO c_vec[j] = int_tabulated( m_vec[0:j], l_vec[0:j] )
     c_vec /= max( c_vec )
     mock_cp[*,i] = interpol( m_vec, c_vec, rnds[*,i] )
  ENDFOR

  If( Keyword_Set( sys_t3 ) ) THEN BEGIN
     rnds_sys = randomu( seed, nSim, n_bispect )
     sys_cp = rnds_sys * 0.0
        For i = 0, n_bispect - 1 DO BEGIN
           
           sys_c = calc_like_dist( (sys_t3.t3phierr)[i], theta_vec=sys_m, mean=(sys_t3.t3phi)[i], /cumulative )
           sys_cp[*,i] = interpol( sys_m, sys_c, rnds_sys[*,i] )
        ENDFOR
     mock_cp += sys_cp
  ENDIF

  ;; Return bispectrum (and this takes care of modulo)
  return, complex( Cos( mock_cp / !radeg ), Sin( mock_cp / !radeg ) )
END
