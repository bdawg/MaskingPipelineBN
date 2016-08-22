;This is an astrometry fitting routine that includes calculating
;the Jacobian for mpfit.pro. This is also used in astrometry_ao_chi2.pro,
;to calculate the chi^2 partial derivatives. 
;NB This differs from astrometry_fit because it doesn't include STEPS
;data, i.e. it only fits for 7 parameters, not 13. A third version of
;this program (RV + AO) should be made...
function astrometry_ao_fit_nojac,  p, data=data,  $
 errors = errors,  rhotheta_ao = rhotheta_ao
 yr =  double(365.256363051)
 nd_ao = n_elements(data.jd_ao)
 pbinary =  p
 rhotheta_ao =  binary_position(pbinary,  data.jd_ao)
 retvect =  [(data.sep_ao-rhotheta_ao[*, 0])/data.sep_err, mod360(data.pa_ao-rhotheta_ao[*, 1])/data.pa_err]
 return,  retvect
end
