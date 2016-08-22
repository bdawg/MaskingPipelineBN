;function to produce a least squares linear fit without a constant term
;cov returns the covariance matrix, but if weights are relative, this must be
;multiplied by MSE, the mean square error. var_yfit should also me multiplied 
;by MSE  

function regress_noc, X,Y,Weights, cov,YFIT, MSE, var_yfit	

On_error,2              ;Return to caller if an error occurs 
SY = SIZE(Y)            ;Get dimensions of x and y.  
SX = SIZE(X)
IF (N_ELEMENTS(Weights) NE SY[1]) OR (SX[0] NE 2) OR (SY[1] NE SX[2]) THEN $
  message, 'Incompatible arrays.'
NTERM = SX[1]           ;# OF TERMS
NPTS = SY[1]            ;# OF OBSERVATIONS
                 
XWY = x # (weights*y)	
WX = fltarr(npts, nterm)
for i = 0,npts-1 do wx(i,*) = x(*,i)*weights(i)
XWX = x # wx
cov = invert(XWX)				;Neter et al pg 402-403
coeff = cov # XWY
yfit = transpose(X) # coeff
if npts ne nterm then MSE = total(weights*(yfit-y)^2)/(npts-nterm)
var_yfit = fltarr(npts)
for i = 0,npts-1 do $
 var_yfit(i) = transpose(X(*,i)) # cov # X(*,i)  ;Neter et al pg 233

return, coeff

end
