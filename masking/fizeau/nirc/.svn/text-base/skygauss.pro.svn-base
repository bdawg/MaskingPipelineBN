pro skygauss,x,a,f,pder
;This is used for Residual Sky Subtraction to match the background...
; X if a vector,A is a parameter vector,F is the function Vector.
;
;For this Gaussian, Let the parameters vector have length 3..
;A=[Area,Displacement,Sigma]
;
;F(X)=Amp * exp (-.5*(x-displacement)^2/(sigma*sigma))
Area=a(0)
disp=a(1)
sigma=a(2)
Amp=Area/(sqrt(2.0*!pi)*sigma)

f=Amp*exp(-.5*(X-disp)*(x-disp)/(sigma*sigma))
pder=fltarr(n_elements(x),3)
for i=0,n_elements(x)-1 do begin
;Area
  pder(i,0)=f(i)/Area
;Disp
  pder(i,1)=f(i)*(x(i)-disp)/(sigma*sigma)
;Sigma
  pder(i,2)=(-1.0*f(i)/sigma)+(f(i)*(x(i)-disp)*(x(i)-disp)/(sigma^3))

endfor
;I hope this doesn't take too long to do..
end
