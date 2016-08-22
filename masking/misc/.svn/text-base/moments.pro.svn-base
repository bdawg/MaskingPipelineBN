;Given a pdf of a distribution(not necessarily normalised), this function 
;returns its first three moments (integral, mean and variance) as a three 
;element vector. If the pdf isn't normalised, each successive moment needs
;to be normalised by dividing my int^(moment)
;The zero point is the position about which moments are to be expanded.

function moments, f, zp=zp

if (keyword_set(zp) eq 0) then zp = 0
n = n_elements(f)
int = total(f)
mn = total((findgen(n) - zp)*f)
mn2 = total((findgen(n) - zp)^2*f)

return, [int,mn,mn2]

end
