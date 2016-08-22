FUNCTION gaussian_1d_lnl, p, data=data, errors=errors

;This is the function that goes along with mmcmc_fit for a (2d gaussian)^2,
;such as for visibility squared.
;
;The parameters are:
; p(0) Amplitude (V_D)
; p(1) alpha (u^2 term)
; p(2) Point Source contribution (V_p) (add a constant)
;
;data is a 4 by n_data -1 array, with
; data[*,0] u coordinates
; data[*,1] src V^2
; data[*,2] cal V^2
;
;errors is a 2 by n_data array, with
; errors[0,*] src V^2 error
; errors[1,*] cal V^2 error

  u = data[*,0]
  model=gaussian_1d_func(u, p)

return, -total(calv2_lnl(model, data[*,1],errors[*,0], data[*,2], errors[*,1]))

END 

