; $Id: //depot/idl/IDL_71/idldir/lib/dist.pro#1 $
;
; Copyright (c) 1982-2009, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;

;+
; NAME:
;	FDIST
;
; PURPOSE:
;	Create a rectangular array in which each element is proportional
;	to its frequency.  This array may be used for a variety
;	of purposes, including frequency-domain filtering and
;	making pretty pictures.
;
; CATEGORY:
;	Signal Processing.
;
; CALLING SEQUENCE:
;	Result = FDIST(N [, M])
;
; INPUTS:
;	N = number of columns in result.
;	M = number of rows in result.  If omitted, N is used to return
;		a square array.
;       xyc = 2-element float array with center coords of required pattern.
;
; OUTPUTS:
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  The computation is done a row at a time.
;
; MODIFICATION HISTORY:
;	Very Old.
; 	SMR, March 27, 1991 - Added the NOZERO keyword to increase efficiency.
;				(Recomended by Wayne Landsman)
;	DMS, July, 1992.  - Added M parameter to make non-square arrays.
;   CT, RSI, March 2000: Changed i^2 to i^2. to avoid overflow.
; PGT hacked up to do distances from fractional pixel location.
; Very crude code just does nested for loops... not fast.
;-
function fdist,n,m,xyc  ;Return a rectangular array in which each pixel = euclidian
		;distance from the origin.
compile_opt idl2

on_error,2              ;Return to caller if an error occurs

n1 = n[0]
m1 = (n_elements(m) le 0) ? n1 : m[0]
 a = FLTARR(n1,m1,/NOZERO)	;Make array

for i=0L, n1-1 do begin	;Row loop
for j=0L, m1-1 do begin	;Row loop
	a[i,j] = sqrt((i-xyc[0])^2 + (j-xyc[1])^2.) ;Euclidian distance
endfor
endfor
return,a
end

