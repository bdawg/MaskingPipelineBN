function conflimits, data, siglist, out

;+
; given a string of data and a list of desired sigma thresholds,
; compute the positive and negative confidence intervals, 
;   spaced about the median (where the median is the returned value)
; 
; will not work right if sample is too small, but it will at least return a warning
; 
; 10/03/07 MCL
;-


if n_params() ne 3 then $
   return, 'function conflimits, data, siglist, out'

message, 'sorting', /info
nn = n_elements(data)
dd = data(sort(data))

; define upper and lower bounds of probability range
plist0 = gauss_pdf(-siglist)
plist1 = gauss_pdf(siglist)

; sanity check
ncheck = 3/(1.0-max(plist1))
if (nn lt ncheck) then begin
   message, '** may not have enough elements in data array for good results **', /info
   message, '     max(sigma) = '+strc(max(siglist)), /info
   message, '     max(probability) = '+strc(max(plist1)), /info
   message, '  desired number of data elements = '+strc(ncheck), /info
   message, '  actual number of data elements = '+strc(nn), /info
endif

message, 'finding', /info
out = fltarr(2, n_elements(siglist))
out(0, *) = dd(round(plist0*nn))
out(1, *) = dd(round(plist1*nn))

return, median(data)
end

