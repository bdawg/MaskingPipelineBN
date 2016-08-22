;This funciton returns the AO (i.e. no STEPS) astrometry chi^2 for inclusion in
;a Monte-Carlo routine. It includes a temperature, T, for simulated-annealing.
function astrometry_ao_chi2,  p,  data=data, chi2_only=chi2_only,$
 errors = errors,  rhotheta_ao = rhotheta_ao,  deriv = deriv,  T = T, $
 pi = pi
 T =  1.
 if (keyword_set(errors)) then T =  errors
 if (keyword_set(pi) eq 0) then begin
  if (data.pi lt 0) then pi = 100.0  $ ;"Typical" parallax in mas
  else pi = data.pi
 endif
 if (p[3] gt 0.999) then addchi2 = 1e9 else addchi2 = 0.
 ;Enforce a prior on eccentricity, to stop circular solutions being
 ;un-naturally preferred.
 addchi2 -= 2*alog(abs(p[3])+1e-3)
 ;Enforce a prior on period and semi-major axis: even in log.
 addchi2 += 2*alog(p[1])
 addchi2 += 2*alog(p[2])
 ;Enforce a prior on inclination: prop proportional to sin(i)
 addchi2 -= 2*alog(abs(sin(p[6]*!pi/180))+1e-4)
 ;Add in the Allen (2007) prior on semi-major axis.
 if (data.aprior[0] ne -1) then $
  addchi2 += (alog10(p[2]/pi) - data.aprior[1])^2/data.aprior[0]^2
 if (data.tprior[0] ne -1) then $
  addchi2 += (alog10(p[1]) - data.tprior[1])^2/data.tprior[0]^2
 ;addchi2 += (alog10(p[2]/pi) - 0.86)^2/0.28^2
 ;addchi2 += (alog10(p[2]/pi) - 0.6)^2/1.0^2
 if (p[2] lt 0) then addchi2 += 1e9
 mod_mtot =  (p[2]/pi)^3/(p[1]/365.25)^2
 addchi2 += (mod_mtot-data.mtot)^2/data.mtot_err^2
 ;addchi2 += ((p[5]-data.epoch)/data.epoch_err)^2
 if (arg_present(deriv)) then begin
  vect =  astrometry_ao_fit(p, jac, data=data, errors = errors,  rhotheta_ao = rhotheta_ao)
  np =  n_elements(p)
  deriv =  fltarr(np)
  for i = 0, np-1 do deriv[i] = 2*total(-vect*jac[*, i])/T
 endif else vect =  astrometry_ao_fit_nojac(p, data=data, errors = errors,  rhotheta_ao = rhotheta_ao)
 chi2_only=total(vect^2)/T
 return, chi2_only + addchi2
end
