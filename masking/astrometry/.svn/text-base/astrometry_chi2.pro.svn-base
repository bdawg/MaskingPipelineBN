;This funciton returns the STEPS+AO astrometry chi^2 for inclusion in
;a Monte-Carlo routine. It includes a temperature, T, for simulated-annealing.
function astrometry_chi2,  p,  data=data, modraoff = modraoff,  moddecoff = moddecoff, $
 errors = errors,  rhotheta_ao = rhotheta_ao,  deriv = deriv,  T = T, chi2_only=chi2_only
 T =  1.
 if (keyword_set(errors)) then T =  errors
 if (p[8] gt 0.999) then addchi2 = 1e9 else addchi2 = 0.
 if (p[12] lt 0) then addchi2 += 1e9
 if (p[7] lt 0) then addchi2 += 1e9
 np =  n_elements(p) < 13
 vect =  astrometry_fit(p, jac, data=data, modraoff = modraoff,  moddecoff = moddecoff, $
  errors = errors,  rhotheta_ao = rhotheta_ao)
if (n_elements(p) gt 14) then begin
  mod_mtot =  (p[12]/(p[4]+p[16]))^3/(p[6]/365.25)^2*(p[12]-p[7])/p[12]
  addchi2 += (mod_mtot-data.mtot)^2/data.mtot_err^2
  addchi2 += (p[16]-data.abspar)^2/data.abspar_sig 
 endif
 ;addchi2 += ((p[5]-data.epoch)/data.epoch_err)^2
 deriv =  fltarr(np)
 for i = 0, np-1 do deriv[i] = 2*total(-vect*jac[*, i])/T
 ;This is special code to include fitting of magnitudes...
 if (n_elements(p) gt 13) then begin
  ;First, add in the astrometry chi^2
  retval = total(vect[0:2*n_elements(data.jd)-1]^2)/T
  ;Now, one epoch at a time, add in the AO chi^2
  for i = 0, n_elements(data.jd_ao)-1 do begin
         if (data.band[i] eq 'J') then model = [rhotheta_ao[i, 0],  rhotheta_ao[i, 1],p[13]] $
    else if (data.band[i] eq 'H') then model = [rhotheta_ao[i, 0],  rhotheta_ao[i, 1],p[14]] $
    else if (data.band[i] eq 'K') then model = [rhotheta_ao[i, 0],  rhotheta_ao[i, 1],p[15]] else stop
    resid = model - [data.sep_ao[i], data.pa_ao[i],  data.crat[i]]
    resid[1] = mod360(resid[1])
    retval += transpose(resid)#invert(data.covar_ao[*, *, i])#resid
  endfor
 endif else retval=total(vect^2)/T
 chi2_only = retval	
 return, chi2_only + addchi2
end
