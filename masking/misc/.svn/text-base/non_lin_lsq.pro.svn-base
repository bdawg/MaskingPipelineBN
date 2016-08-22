;+
; NAME:
;       Non_Lin_Lsq
; PURPOSE:
;       Non-linear least squares fit of a function of an
;       arbitrary number of parameters, to 1-D data set.
;       Function may be any non-linear function where
;       the partial derivatives are known or can be approximated.
; CALLING:
;       Non_Lin_Lsq, Xi, Yd, Parms, sigmas, FUNC_NAME="name of function"
; INPUTS:
;       Xi = vector of independent variables.
;       Yd = vector of dependent variable to be fit with func_name( Xi ), 
;                                                 same length as Xi.
;       Parms = vector of nterms length containing the initial estimate
;              for each parameter.  If Parms is double precision, calculations
;              are performed in double precision, otherwise in single prec.
;              The initial guess of the parameter values should be as close to
;              the actual values as possible or the solution may not converge.
; KEYWORDS:
;       FUNC_NAME = function name (string)
;              Calling mechanism should be:  F = func_name( Xi, Parms, Pderiv )
;         where:
;              F = vector of NPOINT values of function.
;              Xi = vector of NPOINT independent variables, input.
;              Parms = vector of NPARM function parameters, input.
;              Pderiv = array, (NPOINT, NPARM), of partial derivatives.
;                     Pderiv(I,J) = derivative of function at ith point with
;                     respect to jth parameter.  Optional output parameter.
;                     Pderiv should not be calculated if parameter is not
;                     supplied in call (Unless you want to waste some time).
;       WEIGHTS = vector of weights, same length as x and y.
;              For equal (Gaussian) weighting w(i) = 1. (this is default),
;              instrumental (Poisson) weighting w(i) = 1./y(i), etc.
;      /INFO causes Chi-Sq. to be printed each iteration,
;              INFO > 1 causes current parameter estimates to also print.
;       MAX_ITER = maximum # of gradient search iterations, default=20.
;       TOLERANCE = ratio of Chi-Sq. change to previous Chi-Sq. at which
;                     to terminate iterations ( default = 1.e-4 = 0.1% ).
; OUTPUTS:
;       Parms = vector of parameters giving best fit to the data.
; OPTIONAL OUTPUT PARAMETERS:
;       sigmas = Vector of standard deviations for parameters Parms.
;       chisq = final Chi-Sqr deviation of fit.
;       Yfit = resulting best fit to data.
;	CoVariance = covariance matrix for fit params.
; PROCEDURE:
;       Copied from "CURFIT", least squares fit to a non-linear
;       function, pages 237-239, Bevington, Data Reduction and Error
;       Analysis for the Physical Sciences.
;       "This method is the Gradient-expansion algorithm which
;       compines the best features of the gradient search with
;       the method of linearizing the fitting function."
; HISTORY:
;       Written, DMS, RSI, September, 1982.
;       Modified, Frank Varosi, NASA/GSFC, 1992, to use call_function.
;-

pro Non_Lin_Lsq, Xi, Yd, Parms, sigmas, chisq, Yfit, CoVariance, WEIGHTS=Wts, $
	FUNC_NAME=func_name, MAX_ITER=Maxit, INFO_PRINT=info, TOLERANCE=tol

	Nparm = N_ELEMENTS( Parms )			;# of params.
	Ndata = N_elements( Yd )
	Nfree = ( Ndata < N_ELEMENTS(Xi) ) - Nparm     ;degrees of freedom

	IF Nfree LE 0 THEN begin
		message,"not enough data points",/INFO
		return
	   endif

	DIAG = INDGEN( Nparm )*(Nparm+1)	;subscripts of diagonal elements

	if N_elements( Maxit ) NE 1 then begin
		Maxit=20
		notify=1
	  endif else notify=0

	if N_elements( tol ) NE 1 then tol = 1.e-4
	if N_elements( func_name ) NE 1 then func_name = "FUNC"
	if N_elements( Wts ) LT Ndata then Wts = replicate( 1, Ndata )
	Parms = 1.*Parms	;make params floating
	Lambda = 0.001		;Initial lambda
	maxtry = 20

	FOR iter = 1,Maxit DO BEGIN	;Big Iteration loop

		YFIT = call_function( func_name, Xi, Parms, Pderiv )
		BETA = (Yd-YFIT)*Wts # Pderiv
		ALPHA = TRANSPOSE( Pderiv ) # (Wts # (FLTARR(Nparm)+1)*Pderiv)
		ntry=0
		chisq1 = TOTAL( Wts*(Yd-YFIT)^2 )/Nfree   ;present chi squared.

; invert modified curvature matrix to find new parameters.

		REPEAT BEGIN
			C = SQRT( ALPHA(DIAG) # ALPHA(DIAG) )
			ARRAY = ALPHA/C
			ARRAY(DIAG) = 1.+Lambda
			Pnew = Parms + INVERT( ARRAY )/C # TRANSPOSE( BETA )
			Yfit = call_function( func_name, Xi, Pnew )
			chisq = TOTAL( Wts * ( Yd - Yfit )^2 )/Nfree
			ntry = ntry+1
			Lambda = (Lambda*10) < 1.e20	;assume fit got worse

		  ENDREP UNTIL (chisq LE chisq1) OR (ntry GT maxtry)

		Lambda = Lambda/100		;decrease Lambda
		Parms=Pnew			;save new parameter estimate.

		if keyword_set( info ) then begin
			PRINT,iter,"   chisqr =",chisq1,chisq
			if (info GT 1) then PRINT," parms:",Parms
		   endif

		if (chisq EQ 0) then goto,DONE
		IF (chisq1-chisq)/chisq1 LE tol then goto,DONE
	  ENDFOR

	if (notify) OR keyword_set(info) then message,"Failed to converge",/INFO

DONE:	sigmas = SQRT( ARRAY(DIAG) / ALPHA(DIAG) )		;return sigma's

	if N_params() GE 7 then CoVariance = ARRAY/ALPHA
END
