function get_b2, ha, dec, cr=cr

; get_b2	- Calculate the coefficient (b2) in front of x in the
;		dependence of the undistorted coordinates (xp,yp) on the 
;		measured ones (x,y): yp(x,y) = b0+b1*y+b2*x
;		See notes from 03/13/05.
;
; Don't rely on CR dependence.

par = dblarr(3)
a = dblarr(20)
a_err = dblarr(20)

par[0] = ha
par[1] = dec
if n_elements(cr) gt 0 then begin
    par[2] = cr
    a[0] =  1.7948e-07
    a[1] = -3.6233e-06
    a[2] = 1.4093e-07
    a[3] = 6.8511e-11
    a[4] = -5.9924e-06
    a[5] = 2.2882e-07
    a[6] = 1.1431e-09
    a[7] = 1.8815e-07
    a[8] = 9.5323e-09
    a[9] = 5.3582e-10
    a[10] = 2.1818e-07
    a[11] = -4.4001e-06
    a[12] = -1.556e-07
    a[13] = -7.2721e-06
    a[14] = -2.7478e-07
    a[15] = -9.7667e-11
    a[16] = 6.3516e-08
    a[17] = 1.7604e-09
    a[18] = 8.1304e-10
    a[19] = -2.9798e-10
    a_err[0] = 3.8531e-09
    a_err[1] = 9.7206e-08
    a_err[2] = 5.2825e-09
    a_err[3] = 1.0388e-11
    a_err[4] = 7.7077e-08
    a_err[5] = 4.1093e-09
    a_err[6] = 2.33e-11
    a_err[7] = 6.6856e-09
    a_err[8] = 3.0204e-10
    a_err[9] = 7.4936e-11
    a_err[10] = 4.6763e-09
    a_err[11] = 1.1797e-07
    a_err[12] = 4.4868e-09
    a_err[13] = 9.3539e-08
    a_err[14] = 3.5839e-09
    a_err[15] = 5.6545e-11
    a_err[16] = 8.1788e-10
    a_err[17] = 8.894e-11
    a_err[18] = 3.8816e-11
    a_err[19] = 2.3088e-11

    b2 = 0
    k = 0
    for i=0,3 do for j=0,3-i do for m=0,3-i-j do begin 
	    b2 += a[k] * par[0]^i * par[1]^j * par[2]^m
	    k += 1
    endfor
endif else begin
; Based on 15 measurements at zenith, CR = 333.5
;    a[0] =  0.0006344063
;    a[1] = -1.11574e-05
;    a[2] = -6.694108e-08
;    a[3] =  5.538685e-10
;    a[4] =  7.526047e-06
;    a[5] =  0
;    a[6] =  0
;    a[7] =  1.752327e-08
;    a[8] =  7.958071e-10
;    a[9] = -2.83421e-10
;    a_err[0] = 1.958502e-06
;    a_err[1] = 5.531384e-08
;    a_err[2] = 2.372019e-09
;    a_err[3] = 2.402528e-11
;    a_err[4] = 4.054402e-08
;    a_err[5] = 0
;    a_err[6] = 0
;    a_err[7] = 6.278026e-10
;    a_err[8] = 1.136532e-11
;    a_err[9] = 6.209058e-12
; Based on 300 measurements at all hour angles, decs and CR = 333.5
    a[0] =  0.0001996558
    a[1] = -1.054056e-05
    a[2] = -7.407272e-08
    a[3] =  5.37468e-10
    a[4] =  7.580522e-06
    a[5] =  0
    a[6] =  0
    a[7] =  1.065442e-08
    a[8] =  8.868426e-10
    a[9] = -2.833509e-10
    a_err[0] = 5.278138e-06
    a_err[1] = 1.457331e-07
    a_err[2] = 6.395554e-09
    a_err[3] = 6.431538e-11
    a_err[4] = 1.042847e-07
    a_err[5] = 0
    a_err[6] = 0
    a_err[7] = 1.655009e-09
    a_err[8] = 3.064197e-11
    a_err[9] = 1.65091e-11
; Based on 350 measurements at all hour angles, decs and cr angles
;    a[0] =  0.00019503
;    a[1] = -1.0937e-05
;    a[2] = -7.1684e-08
;    a[3] =  5.7651e-10
;    a[4] =  7.6181e-06
;    a[5] =  4.7518e-09
;    a[6] = -7.7406e-11
;    a[7] =  1.6036e-08
;    a[8] =  8.2059e-10
;    a[9] = -2.9566e-10
;    a_err[0] =  6.275e-06
;    a_err[1] = 1.8811e-07
;    a_err[2] = 7.5504e-09
;    a_err[3] = 7.6771e-11
;    a_err[4] =   1.45e-07
;    a_err[5] = 5.4994e-09
;    a_err[6] = 5.6795e-11
;    a_err[7] = 2.2533e-09
;    a_err[8] = 3.8922e-11
;    a_err[9] = 2.3092e-11

    b2 = 0
    k = 0
    for i=0,3 do for j=0,3-i do begin 
	    b2 += a[k] * par[0]^i * par[1]^j
	    k += 1
    endfor
endelse

return, b2

end
