;+
; NAME:
;	gaussian
;
; PURPOSE:
;	Compute the 1-D Gaussian function at an array of points.
;
; CALLING:
;	y = gaussian( xi, parms, pderiv )
;
; INPUTS:
;	xi = array, independent variable of Gaussian function.
;
;	parms = parameters of Gaussian, 2 or 3 element array:
;		parms(0) = maximum value (factor) of Gaussian,
;		parms(1) = mean value (center) of Gaussian,
;		parms(2) = standard deviation (sigma) of Gaussian.
;		parms(3) = optional, constant offset added to Gaussian.
;		(if parms has only 2 elements then sigma taken from common).
;
; OUTPUT:
;	pderiv = optional output of partial derivatives,
;		computed only if parameter is present in call.
;
;		pderiv(*,i) = partial derivative at all xi absisca values
;		with respect to parms(i), i=0,1,2.
;
;	Function returns array of Gaussian evaluated at xi.
;
; COMMON BLOCKS:
;	common gaussian, sigma
; HISTORY:
;	Written: Frank Varosi NASA/GSFC 1992.
;	F.V. 1994, added optional fourth parameter = constant offset.
;-

function gaussian, xi, parms, pderiv

  common gaussian, sigma

	Nparmg = N_elements( parms )
	parms = float( parms )
	if (Nparmg GE 3) then sigma = parms(2)

	z = ( xi - parms(1) )/sigma
	zz = z*z
	gauss = fltarr( N_elements( zz ) )
	w = where( zz LT 172, nw )
	if (nw GT 0) then gauss(w) = exp( -zz(w) / 2 )

	if N_params() GE 3 then begin

		pderiv = fltarr( N_elements( xi ), Nparmg )
		fsig = parms(0) / sigma

		pderiv(0,0) = gauss
		pderiv(0,1) = gauss * z * fsig

		if (Nparmg GE 3) then  pderiv(0,2) = gauss * zz * fsig
		if (Nparmg GE 4) then  pderiv(0,3) = 1
	   endif

	if (Nparmg GE 4) then return, parms(0) * gauss + parms(3) $
			else return, parms(0) * gauss
end
