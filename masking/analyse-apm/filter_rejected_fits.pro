function filter_rejected_fits, start, diff, edir, eyeballer_f, prefix=prefix

  If( eyeballer_f EQ '' ) THEN return, indgen(diff) + start

  If( NOT keyword_set( prefix ) ) THEN prefix = 'ph'
  
  e_f = edir + eyeballer_f
  ; Restores variable called 'rejected_fits'
  restore, e_f

  initarray=1
  okfits = -1 
  For i = FIX( start ), FIX( start ) + FIX( diff ) - 1 DO BEGIN
     fn1 = prefix + String( i, format='(I4.4)' ) + ".fits"
     fn2 = prefix + String( i, format='(I4.4)' ) + ".fits.gz"

     w1 = where( rejected_fits EQ fn1, c1 )
     w2 = where( rejected_fits EQ fn2, c2 )
  
     If( c1 + c2 NE 0 ) THEN CONTINUE
     
     If( initarray EQ 1 ) THEN BEGIN
        okfits = [ i ]
        initarray = 0
     ENDIF ELSE okfits = [ okfits, i ]
  ENDFOR

Return, okfits

END


