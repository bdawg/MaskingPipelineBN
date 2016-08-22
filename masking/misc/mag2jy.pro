; Program to convert input magnitude into 
; Jansky 
;
; Uses photometric calibration of Johnson (1996 ARA&A 4 193)
;
; Input  : mag_vector - vector of mags to be converted to JY
;                    **** MUST BE SORTED to _INCREASING_ Lambda
;        ; wavelengths
; Returns: vector of fluxes, in Jansky
; 
; Still very rough (error trapping etc not done). PGT Sep00

function mag2jy,wavelengths,mag_vector

 
;bands=  U, B, V, R, I, J, H, K, L, M, N
lambda=[0.36, 0.44, 0.55, 0.70, 0.90, 1.25, 1.65, 2.2, 3.4, 5.0, 10.2]
jansky0=[1880., 4440, 3810, 2880, 2240, 1771, 1062, 629, 312, 180, 43]

; NEW set of standards from VEGA (Allen, astrophysical quantities)
; _EXCEPT_ ones indicated in brackets taken from above
; bands= (U),(B),V,(R),(I),J,H,Ks,K,L,L',M,8.7,N,11.7,Q

lambda= [0.36,0.44,0.5556,0.70,0.90,1.215,1.654,2.157,2.179,3.547, $
        3.761,4.769,8.756,10.472,11.653,20.130]
jansky0=[1880,4440,3540,  2880,2240,1630, 1050, 667,  655,  276,   $
        248,  160,  50.0, 35.2,  28.6,  9.70  ]

if(max(wavelengths) gt max(lambda) or min(wavelengths) lt min(lambda)) $
  then print,'### ERROR - Lambda out of interpolation range ### mag2jy(lambdas,mags)'

jansky=spline(lambda,jansky0,wavelengths)


; convert magnitudes
star_flx  = 10^(mag_vector/(-2.5)) * jansky

return,star_flx

end
