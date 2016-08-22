function sky_noise, image, sky_Level, BOX_WIDTH=box_width,  $
					LOCAL_VARIANCE=imvar, $
					HISTO_VARIANCE=imvh,  $
					VALUE_VARIANCE=imvv,  $
					IMAGE_HISTOGRAM=imh,  $
					IMAGE_VALUES=imv,     $
					SIGMA_SET=sigma_set,  $
					debug=debug
;+
; NAME:
;	sky_noise
;
; PURPOSE:
;	Estimate the standard deviation of Gaussian noise by finding the
;	most probable value of standard deviations for moving box of pixels
;	(where the histogram of local standard deviations attains maximum).
;	Optionally, estimate the sky level (average background) by
;	fitting a Gaussian to the histogram of the image.
;
; CALLING:
;	sigma = sky_noise( image, sky_Level, BOX_WIDTH=5 )
;
; INPUTS:
;	image = 2D array of data.
;
; OUTPUTS:
;	sky_Level = the average value of the background in image (optional).
;
;	Function returns the most probable standard deviation of noise in image.
;
; KEYWORDS:
;	BOX_WIDTH = the width of box around each image pixel
;		in which to compute local variances, default = 3.
;
;	SIGMA_SET = if specficied, this value is used as st.dev. of noise
;		and routine skips to determine the sky_Level.
;
; KEYWORD OPTIONAL OUTPUTS:
;
;	LOCAL_VARIANCE = the image of local variances of image (before sqrt)
;	HISTO_VARIANCE = histogram local variance image
;	VALUE_VARIANCE = values corresponding to local variance histogram
;	IMAGE_HISTOGRAM = histogram of image
;	IMAGE_VALUES = values  corresponding to image histogram
;
; EXTERNAL CALLS:
;	function histo
;	pro non_Lin_Lsq
;
; SYSTEM VARIABLES:
;	if !DEBUG GT 3 then histogram plots are produced in temp. window.
; COMMON BLOCKS:
;	common gaussian, sigma_par	;parameter passed to function gaussian.
;	common sky_noise, windeb	;window for debug plots.
; PROCEDURE:
;	Find where the histogram of local standard deviations attains maximum,
;	and iterate changing binsize until it is satisfactorily representing
;	the most probable st.dev., then home in on more exact value
;	by nonlinear Lsq fit to the histogram of local st.devs.
;	Then estimate the sky level (average background) by fitting (the mean)
;	a Gaussian with sigma = most probable st.dev. to the histogram of image.
;
; HISTORY:
;	Written, Frank Varosi NASA/GSFC 1992.
;		 John Monnier ISI 1996  replaced !DEBUG system variable with
;					debug keyword
;-
  common gaussian, sigma_par	;parameter passed to function gaussian.
  common sky_noise, windeb
  if (keyword_set(debug) eq 0) then debug=0 else debug=1
  

	if keyword_set( sigma_set ) then begin
		sigma_noise = sigma_set
		goto,SKY
	   endif
	npix = N_elements( image )
	nbin = fix( sqrt( npix ) ) < 10000
	if N_elements( box_width ) NE 1 then box_width=3

	imvh = histo( image, imvv, NBIN=nbin, BOX_STDEV=box_width, $
							LOCAL_VARIANCE=imvar )
	imvh(0)=0
	maxh = max( imvh, imax )
	imok = 9
	npf = 0.01
	if DEBUG eq 1 then get_window, windeb, XS=500,YS=400
	
	if (maxh LT npix*npf) OR (imax LT imok) then begin

		im_stdev = sqrt( imvar(where( imvar GT 0 )) )
		npix = N_elements( im_stdev )
		nbinp = 0
		iter=0

		while ( (maxh LT npix*npf) OR (imax LT imok) ) AND $
				     (iter LT 20) AND (nbinp NE nbin) do begin
			iter = iter+1
			bfac = 1 + 3/aLog( nbin>2 )
			nbinp = nbin
			if (maxh LT npix*npf) then nbin = nbin/bfac  $
			   else if (imax LT imok) then nbin = bfac*nbin
			nbin = fix( nbin < npix ) > 20
			imvh = histo( im_stdev, imvv, NBIN=nbin )
			imvh(0)=0
			maxh = max( imvh, imax )

			if DEBUG eq 1 then begin
	  			np = ( 3*imax ) < (N_elements( imvv )-1)
				plot,imvv(0:np),imvh(0:np),PSYM=10
				print,iter,imax,maxh,nbin
			   endif
		  endwhile
	   endif

	if dEBUG eq 1 then begin
	  	np = ( 3*imax ) < (N_elements( imvv )-1)
		plot,imvv(0:np),imvh(0:np),PSYM=10,XTIT="St.Devs.",YTIT="Freq."
		print,imax,maxh,nbin
		stop
	   endif
	if (imax GE 2) then begin
	   pg = [ maxh, imvv(imax), (imvv(imax)-imvv(0))/2/sqrt(2*aLog(2)) ]
	   np = ( 3*imax ) < (N_elements( imvv )-1)   ;fit to variance histo

	  np= np < (N_elements( imvh ) -1)

		non_Lin_Lsq, imvv(0:np), imvh(0:np), pg, pgsig, FUNC="gaussian"
		if (pg(1) GT 0) then  sigma_noise = pg(1)  $
				else  sigma_noise = imvv(imax)
	 endif else begin

		sigma_noise = 0
		sky_Level = 0
		return, sigma_noise
	  endelse

	if N_params() LT 2 then return, sigma_noise
SKY:
	sigma_par = sigma_noise	    ;fix sigma parameter of gaussian function.
	imh = histo( image, imv, BINSIZE=(sigma_par/7), NBIN=10000 )

	imh(0)=0				;ignore background pixels
	pg = [ max( imh, imax ), imv(imax) ]	;guess factor and mean,
	np = ( 2*imax ) < (N_elements( imv )-1)	; fit to histogram around Max.
	non_Lin_Lsq, imv(0:np), imh(0:np), pg, pgsig, FUNC="gaussian"
	sky_Level = pg(1)	;estimate of gaussian mean (sky background).

return, sigma_noise
ENd
