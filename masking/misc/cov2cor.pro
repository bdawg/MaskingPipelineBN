;This function returns a correlation matrix when given a covariance
;matrix
;sig returns the standard deviations of the parameters (the sqrt of
;the diagonal)

function cov2cor, cov, sig=sig

n = (size(cov))(1)
if (n ne (size(cov))(2)) then begin
 print, 'Input matrix must be square.'
 stop
endif

x = indgen(n)
sig = sqrt(cov[x,x])
cor = cov
for i = 0,n-1 do for j = 0,n-1 do cor[i,j] = cov[i,j]/sig[i]/sig[j]

return, cor

end
