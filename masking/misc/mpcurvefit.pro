;+
; NAME:
;   MPCURVEFIT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Perform Levenberg-Marquardt least-squares fit (replaces CURVEFIT)
;
; MAJOR TOPICS:
;   Curve and Surface Fitting
;
; CALLING SEQUENCE:
;   YFIT = MPCURVEFIT(X, Y, WEIGHTS, P, [SIGMA,] FUNCTION_NAME=FUNC, 
;                     FUNCTARGS=functargs, ITMAX=itmax, PARINFO=parinfo,
;                     FTOL=ftol, XTOL=xtol, GTOL=gtol, 
;                     ITERPROC=iterproc, ITERARGS=iterargs,
;                     NPRINT=nprint, QUIET=quiet, NOCOVAR=nocovar, 
;                     NFEV=nfev, ITER=iter, ERRMSG=errmsg,
;                     CHISQ=chisq, COVAR=covar, STATUS=status)
;
; DESCRIPTION:
;
;  MPCURVEFIT fits a user-supplied model -- in the form of an IDL
;  function -- to a set of user-supplied data.  MPCURVEFIT calls
;  MPFIT, the MINPACK-1 least-squares minimizer, to do the main
;  work.
;
;  Given the data and their uncertainties, MPCURVEFIT finds the best
;  set of model parameters which match the data (in a least-squares
;  sense) and returns them in the parameter P.  
;
;  MPCURVEFIT returns the best fit function.
;  
;  The user must supply the following items:
;   - An array of independent variable values ("X").
;   - An array of "measured" *dependent* variable values ("Y").
;   - An array of weighting values ("WEIGHTS").
;   - The name of an IDL function which computes Y given X ("FUNC").
;   - Starting guesses for all of the parameters ("P").
;
;  There are very few restrictions placed on X, Y or FUNCT.  Simply
;  put, FUNCT must map the "X" values into "Y" values given the
;  model parameters.  The "X" values may represent any independent
;  variable (not just Cartesian X), and indeed may be multidimensional
;  themselves.  For example, in the application of image fitting, X
;  may be a 2xN array of image positions.
;
;  MPCURVEFIT carefully avoids passing large arrays where possible to
;  improve performance.
;
;  See below for an example of usage.
;   
; INPUTS:
;   FUNCT - a string variable containing the name of an IDL function.
;             This function computes the "model" Y values given the
;             X values and model parameters.  It should be declared in
;             the following way:
;
;             PRO FUNCT, X, parms, YMOD, dparms
;               ; X are the independent variable values
;               ; parms are the parameter values
;               ; dparms - analytical derivative matrix - NOT REQUIRED
;               YMOD   = ...
;               IF N_PARAMS() EQ 3 THEN $
;                 DPARMS = ...  ;; compute derivatives - NOT REQUIRED
;             END
;
;             The returned array YMOD should be of the same size and
;             dimensions as the "measured" Y values.  
;
;             See the discussion of AUTODERIVATIVE and FUNCT in
;             MPFIT.PRO if you wish to compute your derivatives for
;             yourself.  AUTODERIVATIVE is accepted by MPCURVEFIT and
;             passed directly to MPFIT.
;
;   X - Array of independent variable values.
;
;   Y - Array of "measured" dependent variable values.  Y should have
;       the same data type as X.  The function FUNCT should map
;       X->Y.
;
;   WEIGHTS - Array of weights to be used in calculating the
;             chi-squared value.  If WEIGHTS is specified then the ERR
;             parameter is ignored.  The chi-squared value is computed
;             as follows:
;
;                CHISQ = TOTAL( (Y-FUNCT(X,P))^2 * ABS(WEIGHTS) )
;
;             Here are common values of WEIGHTS:
;
;                1D/ERR^2 - Normal weighting (ERR is the measurement error)
;                1D/Y     - Poisson weighting (counting statistics)
;                1D       - Unweighted
;
;   P - An array of starting values for each of the parameters of the
;       model.  The number of parameters should be fewer than the
;       number of measurements.  Also, the parameters should have the
;       same data type as the measurements (double is preferred).
;
;       Upon successful completion the new parameter values are
;       returned in P.
;
;       If both START_PARAMS and PARINFO are passed, then the starting
;       *value* is taken from START_PARAMS, but the *constraints* are
;       taken from PARINFO.
; 
;   SIGMA - The formal 1-sigma errors in each parameter, computed from
;           the covariance matrix.  If a parameter is held fixed, or
;           if it touches a boundary, then the error is reported as
;           zero.
;
; RETURNS:
;
;   Returns the array containing the best-fitting function.
;
; KEYWORD PARAMETERS:
;
;   CHISQ - the value of the summed squared residuals for the
;           returned parameter values.
;
;   COVAR - the covariance matrix for the set of parameters returned
;           by MPFIT.  The matrix is NxN where N is the number of
;           parameters.  The square root of the diagonal elements
;           gives the formal 1-sigma statistical errors on the
;           parameters IF errors were treated "properly" in MYFUNC.
;           Parameter errors are also returned in PERROR.
;
;           To compute the correlation matrix, PCOR, use this:
;           IDL> PCOR = COV * 0
;           IDL> FOR i = 0, n-1 DO FOR j = 0, n-1 DO $
;                PCOR(i,j) = COV(i,j)/sqrt(COV(i,i)*COV(j,j))
;
;           If NOCOVAR is set or MPFIT terminated abnormally, then
;           COVAR is set to a scalar with value !VALUES.D_NAN.
;
;   ERRMSG - a string error or warning message is returned.
;
;   FTOL - a nonnegative input variable. Termination occurs when both
;          the actual and predicted relative reductions in the sum of
;          squares are at most FTOL (and STATUS is accordingly set to
;          1 or 3).  Therefore, FTOL measures the relative error
;          desired in the sum of squares.  Default: 1D-10
;
;   FUNCTARGS - A structure which contains the parameters to be passed
;               to the user-supplied function specified by FUNCT via
;               the _EXTRA mechanism.  This is the way you can pass
;               additional data to your user-supplied function without
;               using common blocks.
;
;               By default, no extra parameters are passed to the
;               user-supplied function.
;
;   GTOL - a nonnegative input variable. Termination occurs when the
;          cosine of the angle between fvec and any column of the
;          jacobian is at most GTOL in absolute value (and STATUS is
;          accordingly set to 4). Therefore, GTOL measures the
;          orthogonality desired between the function vector and the
;          columns of the jacobian.  Default: 1D-10
;
;   ITER - the number of iterations completed.
;
;   ITERARGS - The keyword arguments to be passed to ITERPROC via the
;              _EXTRA mechanism.  This should be a structure, and is
;              similar in operation to FUNCTARGS.
;              Default: no arguments are passed.
;
;   ITERPROC - The name of a procedure to be called upon each NPRINT
;              iteration of the MPFIT routine.  It should be declared
;              in the following way:
;
;              PRO ITERPROC, FUNCT, p, iter, fnorm, FUNCTARGS=fcnargs, $
;                PARINFO=parinfo, QUIET=quiet, ...
;                ; perform custom iteration update
;              END
;         
;              ITERPROC must either accept all three keyword
;              parameters (FUNCTARGS, PARINFO and QUIET), or at least
;              accept them via the _EXTRA keyword.
;          
;              FUNCT is the user-supplied function to be minimized,
;              P is the current set of model parameters, ITER is the
;              iteration number, and FUNCTARGS are the arguments to be
;              passed to FUNCT.  FNORM should be the
;              chi-squared value.  QUIET is set when no textual output
;              should be printed.  See below for documentation of
;              PARINFO.
;
;              In implementation, ITERPROC can perform updates to the
;              terminal or graphical user interface, to provide
;              feedback while the fit proceeds.  If the fit is to be
;              stopped for any reason, then ITERPROC should set the
;              system variable !ERR to a negative value.  In
;              principle, ITERPROC should probably not modify the
;              parameter values, because it may interfere with the
;              algorithm's stability.  In practice it is allowed.
;
;              Default: an internal routine is used to print the
;                       parameter values.
;
;   ITMAX - The maximum number of iterations to perform.  If the
;             number is exceeded, then the STATUS value is set to 5
;             and MPFIT returns.
;             Default: 200 iterations
;
;   NFEV - the number of FUNCT function evaluations performed.
;
;   NOCOVAR - set this keyword to prevent the calculation of the
;             covariance matrix before returning (see COVAR)
;
;   NPRINT - The frequency with which ITERPROC is called.  A value of
;            1 indicates that ITERPROC is called with every iteration,
;            while 2 indicates every other iteration, etc.  Note that
;            several Levenberg-Marquardt attempts can be made in a
;            single iteration.
;            Default value: 1
;
;   PARINFO - Provides a mechanism for more sophisticated constraints
;             to be placed on parameter values.  When PARINFO is not
;             passed, then it is assumed that all parameters are free
;             and unconstrained.  In no case are values in PARINFO
;             modified during a call to MPFIT.
;
;             PARINFO should be an array of structures, one for each
;             parameter.  Each parameter is associated with one
;             element of the array, in numerical order.  The structure
;             can have the following entries (none are required):
;
;               - VALUE - the starting parameter value (but see
;                         START_PARAMS above).
;
;               - FIXED - a boolean value, whether the parameter is to 
;                         be held fixed or not.  Fixed parameters are
;                         not varied by MPFIT, but are passed on to 
;                         FUNCT for evaluation.
;
;               - LIMITED - a two-element boolean array.  If the
;                 first/second element is set, then the parameter is
;                 bounded on the lower/upper side.  A parameter can be
;                 bounded on both sides.  Both LIMITED and LIMITS must
;                 be given together.
;
;               - LIMITS - a two-element float or double array.  Gives
;                 the parameter limits on the lower and upper sides,
;                 respectively.  Zero, one or two of these values can
;                 be set, depending on the value of LIMITED.  Both 
;                 LIMITED and LIMITS must be given together.
;
;               - STEP - the step size to be used in calculating the
;                 numerical derivatives.  If set to zero, then the
;                 step size is computed automatically.
;
;               - TIED - a string expression which "ties" the
;                 parameter to other free or fixed parameters.  Any
;                 expression involving constants and the parameter
;                 array P are permitted.  Example: if parameter 2 is
;                 always to be twice parameter 1 then use the
;                 following: parinfo(2).tied = '2 * P(1)'.  Since they
;                 are totally constrained, tied parameters are
;                 considered to be fixed; no errors are computed for
;                 them.
; 
;             Other tag values can also be given in the structure, but
;             they are ignored.
;
;             Example:
;             parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
;                                  limits:[0.D,0]}, 5)
;             parinfo(0).fixed = 1
;             parinfo(4).limited(0) = 1
;             parinfo(4).limits(0)  = 50.D
;             parinfo(*).value = [5.7D, 2.2, 500., 1.5, 2000.]
;
;             A total of 5 parameters, with starting values of 5.7,
;             2.2, 500, 1.5, and 2000 are given.  The first parameter
;             is fixed at a value of 5.7, and the last parameter is
;             constrained to be above 50.
;
;             Default value:  all parameters are free and unconstrained.
;
;   QUIET - set this keyword when no textual output should be printed
;           by MPFIT
;
;   STATUS - an integer status code is returned.  All values other
;            than zero can represent success.  It can have one of the
;            following values:
;
;	   0  improper input parameters.
;         
;	   1  both actual and predicted relative reductions
;	      in the sum of squares are at most FTOL.
;         
;	   2  relative error between two consecutive iterates
;	      is at most XTOL
;         
;	   3  conditions for STATUS = 1 and STATUS = 2 both hold.
;         
;	   4  the cosine of the angle between fvec and any
;	      column of the jacobian is at most GTOL in
;	      absolute value.
;         
;	   5  the maximum number of iterations has been reached
;         
;	   6  FTOL is too small. no further reduction in
;	      the sum of squares is possible.
;         
;	   7  XTOL is too small. no further improvement in
;	      the approximate solution x is possible.
;         
;	   8  GTOL is too small. fvec is orthogonal to the
;	      columns of the jacobian to machine precision.
;
;   XTOL - a nonnegative input variable. Termination occurs when the
;          relative error between two consecutive iterates is at most
;          XTOL (and STATUS is accordingly set to 2 or 3).  Therefore,
;          XTOL measures the relative error desired in the approximate
;          solution.  Default: 1D-10
;
;
; EXAMPLE:
;
;   ; First, generate some synthetic data
;   npts = 200
;   x  = dindgen(npts) * 0.1 - 10.                  ; Independent variable 
;   yi = gauss1(x, [2.2D, 1.4, 3000.])              ; "Ideal" Y variable
;   y  = yi + randomn(seed, npts) * sqrt(1000. + yi); Measured, w/ noise
;   sy = sqrt(1000.D + y)                           ; Poisson errors
;
;   ; Now fit a Gaussian to see how well we can recover
;   p0 = [1.D, 1., 1000.]                           ; Initial guess
;   yfit = mpcurvefit(x, y, 1/sy^2, p0, $
;                   FUNCTION_NAME='GAUSS1P')        ; Fit a function
;   print, p
;
;   Generates a synthetic data set with a Gaussian peak, and Poisson
;   statistical uncertainty.  Then the same function is fitted to the
;   data to see how close we can get.  GAUSS1 and GAUSS1P are
;   available from the same web page.
;
; REFERENCES:
;
;   MINPACK-1, Jorge More', available from netlib (www.netlib.org).
;   "Optimization Software Guide," Jorge More' and Stephen Wright, 
;     SIAM, *Frontiers in Applied Mathematics*, Number 14.
;
; MODIFICATION HISTORY:
;   Translated from MPFITFUN, 25 Sep 1999, CM
;   Alphabetized documented keywords, 02 Oct 1999, CM
;
;-
; Copyright (C) 1997-1999, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes, and to modify and use for personal or
; internal use, is granted.  All other rights are reserved.
;-

FORWARD_FUNCTION mpcurvefit_eval, mpcurvefit, mpfit

; This is the call-back function for MPFIT.  It evaluates the
; function, subtracts the data, and returns the residuals.
function mpcurvefit_eval, p, dp, _EXTRA=extra

  common mpcurvefit_common, fcn, x, y, wts, f, fcnargs

  ;; The function is evaluated here.  There are four choices,
  ;; depending on whether (a) FUNCTARGS was passed to MPCURVEFIT, which
  ;; is passed to this function as "hf"; or (b) the derivative
  ;; parameter "dp" is passed, meaning that derivatives should be
  ;; calculated analytically by the function itself.
  if n_elements(fcnargs) GT 0 then begin
      if n_params() GT 1 then call_procedure, fcn, x, p, f, dp,_EXTRA=fcnargs $
      else                    call_procedure, fcn, x, p, f,    _EXTRA=fcnargs
  endif else begin
      if n_params() GT 1 then call_procedure, fcn, x, p, f, dp $
      else                    call_procedure, fcn, x, p, f
  endelse

  ;; Compute the deviates, applying the weights
  result = (y-f)*wts
      
  ;; Make sure the returned result is one-dimensional.
  result = reform(result, n_elements(result), /overwrite)
  return, result
  
end

function mpcurvefit, x, y, wts, p, perror, FUNCTARGS=fa, $
                     STATUS=status, QUIET=quiet, $
                     chisq=bestnorm, function_name=fcn, $
                     iter=iter, itmax=maxiter, parinfo=parinfo, $
                     noderivative=noderivative, tol=tol, $
                     covar=covar, errmsg=errmsg, _EXTRA=extra

  status = 0L
  if n_params() EQ 0 then begin
      message, "USAGE: YFIT = MPCURVEFIT(X, Y, WTS, P, DP)", /info
      return, !values.d_nan
  endif
  if n_elements(fcn) EQ 0 then fcn = 'funct'
  if n_elements(noderivative) EQ 0 then noderivative = 0

  common mpcurvefit_common, fc, xc, yc, wc, mc, ac
  fc = fcn & xc = x & yc = y & wc = sqrt(abs(wts)) & mc = 0L
  ac = 0 & dummy = size(temporary(ac))
  if n_elements(fa) GT 0 then ac = fa

  result = mpfit('mpcurvefit_eval', p, maxiter=maxiter, $
                 autoderivative=(1-noderivative), $
                 parinfo=parinfo, STATUS=status, nfev=nfev, BESTNORM=bestnorm,$
                 covar=covar, perror=perror, niter=iter, $
                 ERRMSG=errmsg, quiet=quiet, _EXTRA=extra)

  ;; Retrieve the fit value
  yfit = temporary(mc)
  ;; Now do some clean-up
  xc = 0 & yc = 0 & wc = 0 & mc = 0 & ac = 0

  if NOT keyword_set(quiet) AND errmsg NE '' then $
    message, errmsg, /info $
  else $
    p = result

  return, yfit
end
