;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%               function  wtd_mean.pro                        %
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; Purpose:  to work out weighted mean of sample
;  Original Author -- J.D. Monnier, Unknown date
;  Revision History --  W.C. Danchi, 13Apr99, added code to prevent
;                                             divide by 0 condition

function wtd_mean,vector,error,wdev
ns= size(vector)
weights = fltarr(ns(1),ns(2))
index = where(error le 0.,count)
if count gt 0 then begin
  index2 = where(error gt 0.)
  weights(index2)=1./(error(index2)^2)
  weights(index)=mean(weights(index2))
endif else begin
  weights = 1./(error^2)
endelse

wmean=total(vector*weights)/total(weights)

wdev=1./sqrt(total(weights))

return,wmean
end
