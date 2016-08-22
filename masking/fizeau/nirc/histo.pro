;+
; NAME:
;	Histo
;
; PURPOSE:
;	Get histogram with specified # of bins (or specified binsize),
;	and also output the array and bin values, useful for plotting.
;	Optionally, compute histogram of local standard deviations of data.
;
; CALLING EXAMPLES:
;
;	hist_data = Histo( data, BinVals, NBIN=70, MIN=.1, MAX=50 )
;
;	hist_stdev = Histo( data, BinVals, BOX_STDEV=3 )
;
;	plot, binvals, Histo( data, binvals )
;
; INPUT:
;	data = array of numbers (vector, image, etc.)
;
; KEYWORDS:
;	NBIN = # of bins desired, default=100.
;	BINSIZE = the desired binsize, overrides NBIN in general,
;		but NBIN can still limit the maximum # of bins.
;
;	MIN = minimum data value to include in histogram.
;	MAX = maximum data value to include in histogram
;	     (defaults are actual min or max of data).
;
;	BOX_STDEV = the width in pixels of moving window
;	            in which to compute local variances,
;	            and then histogram of standard deviations is returned.
;
;	LOCAL_VARIANCE = optional output, the array of local variances
;	                 of data (before sqrt), if BOX_STDEV is specified.
;
; OUTPUTS:
;	BinVals = the data values at center of each bin (x-axis for plotting).
;
;	Function returns the histogram of data.
;	If BOX_STDEV is specified then
;	the histogram of local standard deviations of the data is returned.
;
; EXTERNAL CALLS:
;	function vartype
; PROCEDURE:
;	Determine binsize (and bins) and call the IDL histogram function.
; HISTORY:
;	Written, Frank Varosi NASA/GSFC 1989.
;	F.V.1991: added BINSIZE option and then NBIN limits the Max # bins.
;	F.V.1992: added option BOX_STDEV = width of moving box in which to
;	   compute Local standard deviations of data and return histogram of it.
;-

function Histo, data, BinVals, NBIN=Nbin, MIN=minv, MAX=maxv, BINSIZE=binsize, $
					  BOX_STDEV=boxdev, LOCAL_VARIANCE=var

	if keyword_set( boxdev ) then begin
		boxdev = ( 2 * (fix( boxdev )/2) + 1 ) > 3
		bw2 = boxdev^2
		fact = float( bw2 )/(bw2-1)
		var = smooth( (data - smooth( data, boxdev ))^2, boxdev ) *fact
		return, Histo( sqrt( var(where( var GT 0 )) ), BinVals,  $
				BIN=binsize, NBIN=Nbin, MIN=minv, MAX=maxv )
	   endif

	if N_elements( data ) LE 1 then begin
		message,"need an array for data",/INFO
		BinVals = [0,1]
		return,BinVals
	   endif

	if N_elements( maxv ) NE 1 then begin
		maxv = max( data, MIN=mind )
		if N_elements( minv ) NE 1 then  minv = mind
	   endif

	if N_elements( minv ) NE 1 then  minv = min( data )
	if (maxv LE minv) then maxv = 1.1 * minv
	range = double( maxv ) - double( minv )

	if N_elements( binsize ) EQ 1 then begin
		if N_elements( Nbin ) EQ 1 then Maxbin = Nbin
		binsize = double( binsize )
		Nbin = ceil( range/binsize )
		if N_elements( Maxbin ) EQ 1 then begin
			if (Nbin GT Maxbin) then begin
				Nbin = Maxbin
				range = Nbin * binsize
				maxv = range + minv
			   endif
		   endif
	  endif else begin
		if N_elements( Nbin ) NE 1 then Nbin = 101
		binsize = range/(Nbin-1)
	   endelse

	BinVals = binsize * findgen( Nbin ) + minv

	if (vartype( data,/CODE ) LE 3) AND $
	   ( (vartype( binsize,/CODE ) GT 3) OR $
	     (vartype( maxv,/CODE ) GT 3) OR $
	     (vartype( minv,/CODE ) GT 3) ) then begin

	     dens = histogram( float( data ), BIN=binsize, MAX=maxv, MIN=minv )

	 endif else dens = histogram( data, BIN=binsize, MAX=maxv, MIN=minv )

	if N_elements( dens ) GT Nbin then return, dens(0:Nbin-1) $
					else return, dens
end
