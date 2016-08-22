;This funciton returns the AO (i.e. no STEPS) astrometry chi^2 for inclusion in
;a Monte-Carlo routine. It includes a temperature, T, for simulated-annealing.
function astrometry_ao_chi2_ecos_esin,  init_p,  data=data, $
 errors = errors,  rhotheta_ao = rhotheta_ao,  deriv = deriv,  T = T,  pi = pi
 T =  1.
 p = init_p
 p[3] = sqrt(init_p[3]^2+init_p[5]^2)
 p[5] = atan(init_p[3],  init_p[5])*180/!pi
 p[1] = 10.^p[1]
 if (keyword_set(errors)) then T =  errors
 if (keyword_set(pi) eq 0) then pi = 100.0 ;"Typical" parallax in mas
 if (p[3] gt 0.999) then addchi2 = 1e9 else addchi2 = 0.
 if (p[2] lt 0) then addchi2 += 1e9
 np =  n_elements(p)
 vect =  astrometry_ao_fit(p, jac, data=data, errors = errors,  rhotheta_ao = rhotheta_ao)
 mod_mtot =  (p[2]/pi)^3/(p[1]/365.25)^2
 addchi2 += (mod_mtot-data.mtot)^2/data.mtot_err^2
 ;addchi2 += ((p[5]-data.epoch)/data.epoch_err)^2
 deriv =  fltarr(np)
 for i = 0, np-1 do deriv[i] = 2*total(-vect*jac[*, i])/T
 return,  total(vect^2)/T + addchi2
end
