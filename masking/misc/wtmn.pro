;A weighted average routine...
;var returns the variance of the mean, not of x
;MSE is a scaling factor, var(x(i)) = MSE/weights(i)
;Weights are 1/sigma^2: this value of sigma is taken as a worst
;case sigma until there are five data points.
;The error in the mean's variance is
;based on a weighted combination of the MSE and the given sigma,
;with equal weighting for three data points.
;
; original                                                MJI  ????
; changed MSE to option;                                  PGT Nov03
; pass errors rather than weights 
; ignore points with sigman le zero

function wtmn, x, sigma, sdev, MSE=MSE,  wts = wts

ix = where(sigma gt 0.)
if (ix[0] eq -1) then begin
 sdev =  0.0
 MSE = 0.0
 wts = 1.
 return,  x[0]
endif
weights=1./sigma[ix]^2
nelt = float(n_elements(x[ix]))
sumwt = total(weights)
mean = total(x[ix]*weights)/sumwt 
e1 = 1./sumwt
MSE = (nelt gt 1)? total(abs(x[ix]-mean)^2*weights)/(nelt-1):1
;e2 = (nelt gt 4)? MSE:(MSE > 1.0) ;Old formula, incorrect for
;covariances between the input measurements.
e2 = MSE > 1.0 ;ie weights define best possible error
e2 = e2/sumwt
var = (2*e1 + (nelt-1.0)*e2)/(nelt+1.0)
sdev=sqrt(var)

;;Now calculated a weights vector, for the purpose of calculating
;;a covariance matrix elsewhere (e.g. see calibrate_v2_cp.pro)
;;We multiply by sqrt(var/e1) because if var is greater than e1, 
;;the errors have been increased because of the scatter between data
;;sets.
wts = weights/sumwt*sqrt(var/e1)

return, mean
end

