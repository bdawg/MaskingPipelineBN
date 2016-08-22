FUNCTION gauss_splodge_fit_fn, narr, p

;This is the function that returns a circular Gaussian
;
; n is the size of the pixel array for the Gaussian
;The parameters are:
; p(0) Amplitude (V_D)
; p(1) Sigma (width of distribution)
; p(2) Point Source contribution (V_p) (add a constant)
; p(3) x location of peak
; p(4) y location of peak

;
arrsz=(size(narr))[1]

u=fdist(arrsz,arrsz,[p[3],p[4]])
return, model=(p(2)+p(0)*exp(-u^2/p(1)^2))^2 * narr

END 

