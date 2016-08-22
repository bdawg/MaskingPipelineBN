; This routine will fit a baseline offset given
; vectors of ra, dec, fringe_offset

pro fit_bline,ra, dec, fringe_offsetdsfu,vis,err,result=result,sigma=sigma,index=index, $
              ptsrc=ptsrc,guess=guess,power=power

if (keyword_set(sfu) eq 0) then begin
  print,' fit_udisk,sfu,vis,err,params=params,sigma=sigma,result=result,index=index,/ptsrc,guess=guess'
  return
endif

if (keyword_set(power) eq 0) then power=1
power = strtrim(string(power),2)
if (keyword_set(index) eq 0) then index=where(sfu ne -1e-6)
x=sfu(index)
y=vis(index)
er=err(index)

if(keyword_set(ptsrc) ne 0) then goto,pointsrc

if (keyword_set(guess) eq 0) then guess=[mas2rad(20.0),1.0]

; Parameters p(0) diameter in radians
;            p(1) intercept at origin

expr='((2*abs(beselj(3.14159*p(0)*X,1))>1e-6) / ((3.14159*p(0)*X)>1e-6))^'+power+' * p(1)'
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,/quiet)

goto,thatsit

pointsrc:
if (keyword_set(guess) eq 0) then guess=[mas2rad(40.0),0.5,0.5]

; Parameters p(0) diameter in radians
;            p(1) intercept at origin
;            p(2) Point Source Component

expr='((2*abs(beselj(3.14159*p(0)*X,1))>1e-6) / ((3.14159*p(0)*X)>1e-6)) * p(1) + p(2)'
result=mpfitexpr(expr, x, y, er, guess, perror=sigma,nprint=100)

thatsit:

end
