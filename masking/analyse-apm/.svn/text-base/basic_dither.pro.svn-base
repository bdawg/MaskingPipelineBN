pro basic_dither, stars, mask, sky, clean_stars, smooth_radius=smooth_radius, percent_use_median=percent_use_median

  If( NOT keyword_set( smooth_radius ) ) THEN smooth_radius = 15.
  If( NOT keyword_set( percent_use_median ) ) THEN percent_use_median = .10

  sz = size( stars, /dimension )
  naxis1 = sz[0]
  naxis2 = sz[1]
  nimages = sz[2]

  If( n_elements( sz ) NE 3 ) THEN stop

  ;; Aside, expand mask to the size of images
  sz = (size(mask, /dimensions ) )
  maxis1 = sz[0]
  maxis2 = sz[1]
  full_mask = 1.0 + fltarr( naxis1, naxis2 )
  full_mask[0:(maxis1-1), 0:(maxis2-1)] = mask

  ;; The remaining flux is 'probably' signal.
  ;; Second, for each image, determine the location of the mask that
  ;; covers the most remaining flux.

  mask_loc = fltarr( nimages, 2 )
  For i = 0, nimages - 1 DO BEGIN
     smooth_image = smooth( stars[ *, *, i ], smooth_radius, /edge_truncate )
     dummy = Max( smooth_image, max_1dloc )
     ; convert to 2D indices
     mask_loc[i,*] = array_indices( smooth_image, max_1dloc )
  ENDFOR
  
  ;; Third, with the mask in place on each frame, calculate the median
  ;; of the non-masked pixels.  This is the true background.
  
  masked_star = fltarr( naxis1, naxis2, nimages )
  For i = 0, nimages - 1 DO BEGIN
     this_mask = shift( full_mask, mask_loc[i,0] - maxis1/2, mask_loc[i,1] - maxis2/2 )
     masked_star[*,*,i] = stars[*,*,i] * this_mask 
  ENDFOR

  sky = fltarr( naxis1, naxis2 )
  clean_stars = fltarr( naxis1, naxis2, nimages )
  For i = 0, naxis1 - 1 DO For j = 0, naxis2 - 1 DO BEGIN
     w = where( Finite( masked_star[i,j,*] ), count )
     If( count GE percent_use_median * nimages ) THEN sky[i,j] = median( (stars[i,j,*])[w] ) ELSE sky[i,j]=!values.f_nan
  ENDFOR

   ; If there are any pixels for which there is no sky data (NaN above),
   ; set it equal to the median sky
  
  w = where( Finite( sky ), complement=needed, goodcount ) 
  med_sky = median( sky[w] )
  If( goodcount NE n_elements( sky ) ) THEN sky[needed] = med_sky

  For i = 0, nimages-1 DO clean_stars[*,*,i] = stars[*,*,i] - sky

END
