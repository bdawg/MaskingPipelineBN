; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    uniform_disk                                %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Based on display_vis.pro (PGT)
;	    Wrote Program.				JDM 12/2/96
;           added power option for v^2                  PGT Nov03
; %%%%% This is a function to return the visibility of a uniform
;       disk as a function of baseline.
;	Note: F is positive definite
;
; To Add: Do partial Derivatives EXplicitly.

pro uniform_disk,x,a,f,pder,power=power  ;pder not supported

; There are only two parameters to this fit.
;      a(0)=The Diameter in inverse units of x.
; Note a negative diameter will return 1/J1.. which is as if 
;  one reversed the roles of calibrator and source.
;      a(1)=The zero visibility intercept.  Generally should be 1.0
;
if (keyword_set(power) eq 0) then power=1

pi=3.1415926

diameter=abs(a(0))
intercept=a(1)
index=where(x eq 0,count)
if (count gt 0) then x(index)=1e-15

f=2*beselj(pi*diameter*x,1) / (pi*diameter*x)

if (a(0) ge 0) then f=abs(intercept*f^power) else f=abs(intercept/f^power)

; AT the MOMENT, no pder calculation.


end
