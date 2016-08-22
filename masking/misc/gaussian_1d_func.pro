FUNCTION gaussian_1d_func, u, p

;This is the function that returns a circular Gaussian
;
;The parameters are:
; p(0) Amplitude (V_D)
; p(1) alpha and beta (u^2 and v^2 term)
; p(2) Point Source contribution (V_p) (add a constant)
;
return, model=(p(2)+p(0)*exp(-p(1)*u^2))^2

END 

