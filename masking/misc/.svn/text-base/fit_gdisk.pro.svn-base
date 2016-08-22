; This routine will fit a Gaussian Disk model to a set of data.
; Peter Tuthill March 6, 1999
; added ability to force origin to equal force0 (probably 1.0)         PGT Nov03
;   Note: pre Nov03 version of mpfitexpr does _not_ work for 1 param fits
;
pro fit_gdisk,sfu,vis,err,result=result,sigma=sigma,index=index, $
              ptsrc=ptsrc,force0=force0,guess=guess,power=power

if (keyword_set(sfu) eq 0) then begin
  print,' fit_gdisk,sfu,vis,err,params=params,sigma=sigma,result=result,index=index,/ptsrc,guess=guess'
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
if(keyword_set(ptsrc) ne 0) then goto,pointsrc

if (keyword_set(guess) eq 0) then guess=[mas2rad(10.0),1.0]
; Parameters p(0) diameter in radians
;            p(1) intercept at origin
expr='exp(-3.57*'+power+'*p(0)^2*x^2) * p(1)'
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,/quiet)
goto,thatsit


pointsrc:
if (keyword_set(guess) eq 0) then guess=[mas2rad(10.0),0.5,0.5]
; Parameters p(0) diameter in radians
;            p(1) intercept at origin
;            p(2) Point Source Component
expr='exp(-3.57*'+power+'*p(0)^2*x^2) * p(1) + p(2)'
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,nprint=100)
goto,thatsit


oneparamfit: ;diam only, no intercept fit
if (keyword_set(guess) eq 0) then guess=[mas2rad(10.0),1.0]
parinfo = replicate({fixed:0},2)
parinfo(1).fixed = 1   ;this fixes P(1) so it can't vary
expr='exp(-3.57*'+power+'*p(0)^2*x^2) * p(1)'
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,/quiet,parinfo=parinfo)
goto,thatsit

thatsit:

end
