function get_a1, ha, dec, cr=cr

; get_a1	- Calculate the coefficient (a1) in front of y in the
;		dependence of the undistorted coordinates (xp,yp) on the 
;		measured ones (x,y): xp(x,y) = a0+a1*y+a2*y^2+a3*x+a5*x^2
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
    a[0] = -2.4674e-07
    a[1] =  5.6345e-06
    a[2] = -2.6232e-07
    a[3] = -3.4738e-11
    a[4] = 4.7395e-06
    a[5] = -1.3549e-07
    a[6] = -1.0035e-09
    a[7] = -3.4557e-07
    a[8] = -1.653e-08
    a[9] = -6.0375e-10
    a[10] = -2.9955e-07
    a[11] = 6.8375e-06
    a[12] = 2.4715e-07
    a[13] = 5.7525e-06
    a[14] = 2.18e-07
    a[15] = 2.7295e-10
    a[16] = -5.0532e-08
    a[17] = -2.5936e-09
    a[18] = -9.9752e-10
    a[19] = 3.6323e-10
    a_err[0] = 3.885e-09
    a_err[1] = 9.8139e-08
    a_err[2] = 5.3265e-09
    a_err[3] = 1.0587e-11
    a_err[4] = 7.7852e-08
    a_err[5] = 4.1089e-09
    a_err[6] = 2.3212e-11
    a_err[7] = 6.8501e-09
    a_err[8] = 3.0671e-10
    a_err[9] = 7.48e-11
    a_err[10] = 4.7131e-09
    a_err[11] = 1.1903e-07
    a_err[12] = 4.5248e-09
    a_err[13] = 9.4462e-08
    a_err[14] = 3.6168e-09
    a_err[15] = 5.6326e-11
    a_err[16] = 8.2537e-10
    a_err[17] = 8.9195e-11
    a_err[18] = 3.8689e-11
    a_err[19] = 2.3062e-11

    a1 = 0
    k = 0
    for i=0,3 do for j=0,3-i do for m=0,3-i-j do begin 
	    a1 += a[k] * par[0]^i * par[1]^j * par[2]^m
	    k += 1
    endfor
endif else begin
; Based on 15 measurements at zenith, CR = 333.5
;    a[0] = -0.0003955622
;    a[1] =  7.237992e-06
;    a[2] =  8.954493e-08
;    a[3] = -5.357639e-10
;    a[4] = -8.125414e-06
;    a[5] =  0
;    a[6] =  0
;    a[7] =  1.745292e-08
;    a[8] = -9.679733e-10
;    a[9] =  3.561776e-10
;    a_err[0] =  2.048394e-06
;    a_err[1] =  5.58073e-08
;    a_err[2] =  2.418647e-09
;    a_err[3] =  2.461147e-11
;    a_err[4] =  4.214094e-08
;    a_err[5] =  0
;    a_err[6] =  0
;    a_err[7] =  6.626249e-10
;    a_err[8] =  1.195152e-11
;    a_err[9] =  6.472467e-12
; Based on 300 measurements at all hour angles, decs and CR = 333.5
    a[0] = -0.0002496325
    a[1] =  7.222127e-06
    a[2] =  9.158649e-08
    a[3] = -5.467939e-10
    a[4] = -8.138857e-06
    a[5] =  0
    a[6] =  0
    a[7] =  2.137211e-08
    a[8] = -1.030719e-09
    a[9] =  3.659832e-10
    a_err[0] =  5.378355e-06
    a_err[1] =  1.469128e-07
    a_err[2] =  6.465733e-09
    a_err[3] =  6.509733e-11
    a_err[4] =  1.060632e-07
    a_err[5] =  0
    a_err[6] =  0
    a_err[7] =  1.687558e-09
    a_err[8] =  3.126475e-11
    a_err[9] =  1.681164e-11
; Based on 350 measurements at all hour angles, decs and CR angles
;    a[0] = -0.00033358
;    a[1] =  7.4561e-06
;    a[2] =  9.6272e-08
;    a[3] = -6.2675e-10
;    a[4] = -8.0115e-06
;    a[5] = -2.1387e-08
;    a[6] =  2.6401e-10
;    a[7] =  1.8619e-08
;    a[8] =  -1.001e-09
;    a[9] =  3.6227e-10
;    a_err[0] =  6.3316e-06
;    a_err[1] =   1.875e-07
;    a_err[2] =  7.5547e-09
;    a_err[3] =  7.6548e-11
;    a_err[4] =  1.4491e-07
;    a_err[5] =  5.4916e-09
;    a_err[6] =  5.6549e-11
;    a_err[7] =  2.2583e-09
;    a_err[8] =  3.8793e-11
;    a_err[9] =  2.3065e-11

    a1 = 0
    k = 0
    for i=0,3 do for j=0,3-i do begin 
	    a1 += a[k] * par[0]^i * par[1]^j
	    k += 1
    endfor
endelse

return, a1

end
