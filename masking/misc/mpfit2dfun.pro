;+
; NAME:
;   MPFIT2DFUN
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Perform Levenberg-Marquardt least-squares fit to a 2-D IDL function
;
; MAJOR TOPICS:
;   Curve and Surface Fitting
;
; CALLING SEQUENCE:
;   parms = MPFIT2DFUN(MYFUNCT, X, Y, Z, ERR, start_parms, ...)
;
; DESCRIPTION:
;
;  MPFIT2DFUN fits a user-supplied model -- in the form of an IDL
;  function -- to a set of user-supplied data.  MPFIT2DFUN calls
;  MPFIT, the MINPACK-1 least-squares minimizer, to do the main
;  work.  MPFIT2DFUN is a specialized version for two-dimensional 
;  data.
;
;  Given the data and their uncertainties, MPFIT2DFUN finds the best set
;  of model parameters which match the data (in a least-squares
;  sense) and returns them in an array.
;  
;  The user must supply the following items:
;   - Two arrays of independent variable values ("X", "Y").
;   - An array of "measured" *dependent* variable values ("Z").
;   - An array of "measured" 1-sigma uncertainty values ("ERR").
;   - The name of an IDL function which computes Z given (X,Y) ("MYFUNCT").
;   - Starting guesses for all of the parameters ("START_PARAMS").
;
;  There are very few restrictions placed on X, Y, Z, or MYFUNCT.
;  Simply put, MYFUNCT must map the (X,Y) values into Z values given
;  the model parameters.  The (X,Y) values are usually the independent
;  X and Y coordinate positions in the two dimensional plane, but need
;  not be.
;
;  MPFIT2DFUN carefully avoids passing large arrays where possible to
;  improve performance.
;
;  See below for an example of usage.
;   
; INPUTS:
;   MYFUNCT - a string variable containing the name of an IDL
;             function.  This function computes the "model" Z values
;             given the X,Y values and model parameters.  It should be
;             declared in the following way:
;
;             FUNCTION MYFUNCT, X, Y, parms, dparms
;               ; X,Y are the independent variable values
;               ; parms are the parameter values
;               ; dparms - analytical derivative matrix - NOT REQUIRED
;               ZMOD   = ...
;               IF N_PARAMS() EQ 3 THEN $
;                 DPARMS = ...  ;; compute derivatives - NOT REQUIRED
;               RETURN, ZMOD
;             END
;
;             The returned array ZMOD should be of the same size and
;             dimensions as the "measured" Z values.  
;
;             See the discussion of AUTODERIVATIVE and MYFUNCT in
;             MPFIT.PRO if you wish to compute your derivatives for
;             yourself.  AUTODERIVATIVE is accepted by MPFIT2DFUN and
;             passed directly to MPFIT.
;
;   X - Array of "X" independent variable values.  These values are
;       passed directly to the fitting function unmodified.  For image
;       fitting applications, this variable should be a
;       two-dimensional *array* describing the X position of every
;       *pixel*.  [ Thus irregular sampling is permitted. ]  If the
;       sampling is regular, and your x- and y-coordinates are labeled
;       by XP and YP, then you can compute these arrays as follows.
;
;           X = XP # (YP*0 + 1)  ;; Uniform X values
;           Y = (XP*0 + 1) # YP  ;; Uniform Y values
;
;   Y - Array of "Y" indendent variable values.  See discussion under
;       X above.  X and Y should have the same dimension and data
;       type.
;
;   Z - Array of "measured" dependent variable values.  Z should have
;       the same data type as X and Y.  The function MYFUNCT should
;       map (X,Y)->Z.
;
;   ERR - Array of "measured" 1-sigma uncertainties.  ERR should have
;         the same data type as Z.  ERR is ignored if the WEIGHTS
;         keyword is specified.
;
;   START_PARAMS - An array of starting values for each of the
;                  parameters of the model.  The number of parameters
;                  should be fewer than the number of measurements.
;                  Also, the parameters should have the same data type
;                  as the measurements (double is preferred).
;
;                  This parameter is optional if the PARINFO keyword
;                  is used (see MPFIT).  The PARINFO keyword provides
;                  a mechanism to fix or constrain individual
;                  parameters.  If both START_PARAMS and PARINFO are
;                  passed, then the starting *value* is taken from
;                  START_PARAMS, but the *constraints* are taken from
;                  PARINFO.
; 
; RETURNS:
;
;   Returns the array of best-fit parameters.
;
; KEYWORD PARAMETERS:
;
;   BESTNORM - the value of the summed squared residuals for the
;              returned parameter values.
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
;               to the user-supplied function specified by MYFUNCT via
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
;   ITERARGS - The keyword arguments to be passed to ITERPROC via the
;              _EXTRA mechanism.  This should be a structure, and is
;              similar in operation to FUNCTARGS.
;              Default: no arguments are passed.
;
;   ITERPROC - The name of a procedure to be called upon each NPRINT
;              iteration of the MPFIT routine.  It should be declared
;              in the following way:
;
;              PRO ITERPROC, MYFUNCT, p, iter, fnorm, FUNCTARGS=fcnargs, $
;                PARINFO=parinfo, QUIET=quiet, ...
;                ; perform custom iteration update
;              END
;         
;              ITERPROC must either accept all three keyword
;              parameters (FUNCTARGS, PARINFO and QUIET), or at least
;              accept them via the _EXTRA keyword.
;          
;              MYFUNCT is the user-supplied function to be minimized,
;              P is the current set of model parameters, ITER is the
;              iteration number, and FUNCTARGS are the arguments to be
;              passed to MYFUNCT.  FNORM should be the
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
;   MAXITER - The maximum number of iterations to perform.  If the
;             number is exceeded, then the STATUS value is set to 5
;             and MPFIT returns.
;             Default: 200 iterations
;
;   NFEV - the number of MYFUNCT function evaluations performed.
;
;   NITER - the number of iterations completed.
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
;                         MYFUNCT for evaluation.
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
;   PERROR - The formal 1-sigma errors in each parameter, computed
;            from the covariance matrix.  If a parameter is held
;            fixed, or if it touches a boundary, then the error is
;            reported as zero.
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
;   WEIGHTS - Array of weights to be used in calculating the
;             chi-squared value.  If WEIGHTS is specified then the ERR
;             parameter is ignored.  The chi-squared value is computed
;             as follows:
;
;                CHISQ = TOTAL( (Z-MYFUNCT(X,Y,P))^2 * ABS(WEIGHTS) )
;
;             Here are common values of WEIGHTS:
;
;                1D/ERR^2 - Normal weighting (ERR is the measurement error)
;                1D/Z     - Poisson weighting (counting statistics)
;                1D       - Unweighted
;
;   XTOL - a nonnegative input variable. Termination occurs when the
;          relative error between two consecutive iterates is at most
;          XTOL (and STATUS is accordingly set to 2 or 3).  Therefore,
;          XTOL measures the relative error desired in the approximate
;          solution.  Default: 1D-10
;
;   YFIT - the best-fit model function, as returned by MYFUNCT.
;
; EXAMPLE:
;
;   p  = [2.2D, -0.7D, 1.4D, 3000.D]
;   x  = (dindgen(200)*0.1 - 10.) # (dblarr(200) + 1)
;   y  = (dblarr(200) + 1) # (dindgen(200)*0.1 - 10.)
;   zi = gauss2(x, y, p)
;   sz = sqrt(zi)
;   z  = zi + randomn(seed, 200, 200) * sz
;
;   p0 = [0D, 0D, 1D, 10D]
;   p = mpfit2dfun('GAUSS2', x, y, z, sz, p0)
;   
;   Generates a synthetic data set with a Gaussian peak, and Poisson
;   statistical uncertainty.  Then the same function (but different
;   starting parameters) is fitted to the data to see how close we can
;   get.
;
;   It is especially worthy to notice that the X and Y values are
;   created as full images, so that a coordinate is attached to each
;   pixel independently.  This is the format that GAUSS2 accepts, and
;   the easiest for you to use in your own functions.
;
; REFERENCES:
;
;   MINPACK-1, Jorge More', available from netlib (www.netlib.org).
;   "Optimization Software Guide," Jorge More' and Stephen Wright, 
;     SIAM, *Frontiers in Applied Mathematics*, Number 14.
;
; MODIFICATION HISTORY:
;   Written, transformed from MPFITFUN, 26 Sep 1999, CM
;   Alphabetized documented keywords, 02 Oct 1999, CM
;   Added example, 02 Oct 1999, CM
;
;-
; Copyright (C) 1997-1999, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes, and to modify and use for personal or
; internal use, is granted.  All other rights are reserved.
;-

FORWARD_FUNCTION mpfit2dfun_eval, mpfit2dfun, mpfit

; This is the call-back function for MPFIT.  It evaluates the
; function, subtracts the data, and returns the residuals.
function mpfit2dfun_eval, p, dp, _EXTRA=extra

  common mpfit2dfun_common, fcn, x, y, z, err, wts, f, fcnargs

  ;; The function is evaluated here.  There are four choices,
  ;; depending on whether (a) FUNCTARGS was passed to MPFIT2DFUN, which
  ;; is passed to this function as "hf"; or (b) the derivative
  ;; parameter "dp" is passed, meaning that derivatives should be
  ;; calculated analytically by the function itself.
  if n_elements(fcnargs) GT 0 then begin
      if n_params() GT 1 then f = call_function(fcn,x,y,p, dp, _EXTRA=fcnargs)$
      else                    f = call_function(fcn,x,y,p,     _EXTRA=fcnargs)
  endif else begin
      if n_params() GT 1 then f = call_function(fcn,x,y,p, dp) $
      else                    f = call_function(fcn,x,y,p)
  endelse

  ;; Compute the deviates, applying either errors or weights
  if n_elements(err) GT 0 then begin
      result = (z-f)/err
  endif else if n_elements(wts) GT 0 then begin
      result = (z-f)*wts
  endif else begin
      result = (z-f)
  endelse
      
  ;; Make sure the returned result is one-dimensional.
  result = reform(result, n_elements(result), /overwrite)
  return, result
  
end

function mpfit2dfun, fcn, x, y, z, err, p, WEIGHTS=wts, FUNCTARGS=fa, $
                   BESTNORM=bestnorm, nfev=nfev, STATUS=status, $
                   parinfo=parinfo, $
                   covar=covar, perror=perror, niter=iter, yfit=yfit, $
                   quiet=quiet, ERRMSG=errmsg, _EXTRA=extra

  status = 0L
  if n_params() EQ 0 then begin
      message, "USAGE: PARMS = MPFIT2DFUN('MYFUNCT', X, Y, ERR, "+ $
        "START_PARAMS, ... )", /info
      return, !values.d_nan
  endif

  ;; Use common block to pass data back and forth
  common mpfit2dfun_common, fc, xc, yc, zc, ec, wc, mc, ac
  fc = fcn & xc = x & yc = y & zc = z & mc = 0L
  ;; These optional parameters must be undefined first
  ac = 0 & dummy = size(temporary(ac))
  ec = 0 & dummy = size(temporary(ec))
  wc = 0 & dummy = size(temporary(wc))

  if n_elements(fa) GT 0 then ac = fa
  if n_elements(wts) GT 0 then begin
      wc = sqrt(abs(wts))
  endif else if n_elements(err) GT 0 then begin
      wh = where(err EQ 0, ct)
      if ct GT 0 then begin
          message, 'ERROR: ERROR value must not be zero.  Use WEIGHTS.', $
            /info
          return, !values.d_nan
      endif
      ec = err
  endif

  result = mpfit('mpfit2dfun_eval', p, $
                 parinfo=parinfo, STATUS=status, nfev=nfev, BESTNORM=bestnorm,$
                 covar=covar, perror=perror, niter=iter, $
                 ERRMSG=errmsg, quiet=quiet, _EXTRA=extra)

  ;; Retrieve the fit value
  yfit = temporary(mc)
  ;; Some cleanup
  xc = 0 & yc = 0 & zc = 0 & wc = 0 & ec = 0 & mc = 0 & ac = 0

  ;; Print error message if there is one.
  if NOT keyword_set(quiet) AND errmsg NE '' then $
    message, errmsg, /info

  return, result
end
