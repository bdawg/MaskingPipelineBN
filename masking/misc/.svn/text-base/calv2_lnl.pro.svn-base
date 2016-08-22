;This is Mike's likelihood function...
;Given a src V^2 + error and a cal V^2 + error 
;this calculates a non-normalised log-likelihood for 
;calibrated V^2

function calv2_lnl, xin, mus,sigs, muc, sigc

x = xin
;Normalise everything to sigc...
sigs = sigs/sigc
mus = mus/sigc
muc = muc/sigc

;Now normalise mus to sigs
mus = mus/sigs
x = x/sigs
;temp = (mus*x+muc)^2/2.0/(1.0+x^2)
;retarr = x

;w1 = where(temp gt 10.0, complement=w2)
;if (w1[0] ne -1) then  retarr[w1] = $
; temp[w1] + alog((mus*x[w1]+muc)/(1+x[w1]^2)^1.5*erf((mus*x[w1] +muc)/sqrt(2.0*(x[w1]^2+1.0)))) 
;if (w2[0] ne -1) then  retarr[w2] = $
; alog(sqrt(2/!pi)*(1+x[w2]^2) + (mus*x[w2]+muc)*erf((mus*x[w2] +muc)/sqrt(2.0*(x[w2]^2+1.0)))*exp(temp[w2])) - 1.5*alog(1+x[w2]^2)

;return, retarr
;brak = (mus*x +muc)/(x^2+1.0)

;return, sqrt(x^2/(1+x^2))*exp(-(mus-muc*x)^2/2/(1+x^2))
;return, sqrt(1.0/(1+x^2)/2.0/!pi)*exp(-(mus^2+muc^2)/2 + (mus*x +muc)^2/2/(1+x^2))
;return, (mus*x +muc)^2/2.0/(1.0+x^2);-1.5*alog(1.0+x^2) +
;alog(mus*x+muc)

return, -alog(2*!pi*sigs) - (muc^2+mus^2)/2.0 + (mus*x +muc)^2/2.0/(1.0+x^2)+ $
  alog(2*exp(-(mus*x +muc)^2/(x^2+1.0)/2)/(x^2+1.0) + sqrt(2.0*!pi)*(mus*x +muc)/(x^2+1)^1.5*erf((mus*x +muc)/sqrt((x^2+1)*2.0)))
;return, alog( 1./(1+x^2)*(sqrt(2.0) + sqrt(!pi)*(mus*x+muc)/sqrt(1+x^2)*erf((mus*x +muc)/sqrt(2.0*(x^2+1)))* $
;              exp((mus*x+muc)^2/2.0/(1.0+x^2))) )
;return, (mus*x+muc)^2/2.0/(1.0+x^2) + alog((mus*x+muc)/(1+x^2)^(1.5)*erf((mus*x +muc)/sqrt(2.0*(x^2+1.0))))

end
