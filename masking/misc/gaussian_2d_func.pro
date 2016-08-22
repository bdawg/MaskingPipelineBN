FUNCTION gaussian_2d_func, p, u, v

;This is the function that goes along with mmcmc_fit for a (2d gaussian)^2,
;such as for visibility squared.
;
;The parameters are:
; p(0) Amplitude (V_D)
; p(1) alpha (u^2 term)
; p(2) beta (v^2 term)
; p(3) gamma (cross uv term)
; p(4) Point Source contribution (V_p) (add a constant)
;
return, model=(p(4)+p(0)*exp(-(p(1)*u^2+p(2)*v^2+p(3)*u*v)))^2

END 

