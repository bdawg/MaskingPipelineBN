; This routine will fit a model of 2 Gaussian Disks to a set of data.
; Peter Tuthill March 6, 1999
; added ability to force origin to equal force0 (probably 1.0)         PGT Nov03
;   Note: pre Nov03 version of mpfitexpr does _not_ work for 1 param fits
;
pro fit_2gdisk,sfu,vis,err,result=result,sigma=sigma,index=index, $
              force0=force0,guess=guess,power=power

if (keyword_set(sfu) eq 0) then begin
  print,' fit_2gdisk,sfu,vis,err,params=params,sigma=sigma,result=result,index=index,guess=guess'
  return
endif

if (keyword_set(power) eq 0) then power=1
power = strtrim(string(power),2)
if (keyword_set(index) eq 0) then index=where(sfu ne -1e-6)
x=sfu(index)
y=vis(index)
er=err(index)

if (keyword_set(force0) eq 0) then force0=0
if(force0 ne 0) then goto,oneparamfit

if (keyword_set(guess) eq 0) then guess=[mas2rad(10.0),0.5,3.0,0.5]
; Parameters p(0) diameter component 1 in radians
;            p(1) flux, component 1 
;            p(2) diameter component 2 in radians
;            p(3) flux, component 2 

expr1='exp(-3.57*p(0)^2*x^2) * p(1)'
expr2='exp(-3.57*p(2)^2*x^2) * p(3)'
expr='('+expr1+' + '+expr2+')^'+power

result=mpfitexpr(expr, x, y, er, guess, perror=sigma,/quiet)
goto,thatsit

oneparamfit: ; 3 free params actually: fwhm1, 2 and flux ratio (total=1)
if (keyword_set(guess) eq 0) then guess=[mas2rad(10.0),0.5,3.0,0.5]
parinfo = replicate({tied:''},4)
parinfo(3).tied = '1.0-P(1)'   ;this fixes P(1) + P(3) = 1.0
expr1='exp(-3.57*p(0)^2*x^2) * p(1)'
expr2='exp(-3.57*p(2)^2*x^2) * p(3)'
expr='('+expr1+' + '+expr2+')^'+power
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,/quiet,parinfo=parinfo)
goto,thatsit

thatsit:

end
