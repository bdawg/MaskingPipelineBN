;This is an astrometry fitting routine that includes the calculating
;the Jacobian for mpfit.pro. This is also used in astrometry_chi2.pro,
;to calculate the chi^2 partial derivatives. 
;TODO: -As this has lots of improvements since older masking/SUSI orbit
; fitting, I should really make a no-STEPS version of the function.
function astrometry_fit,  p,  jac, data=data, modraoff = modraoff,  moddecoff = moddecoff, $
 errors = errors,  rhotheta_ao = rhotheta_ao,  bin_ra = bin_ra,  bin_dec = bin_dec
 yr =  double(365.256363051)
 nd = n_elements(data.jd)
 nd_ao = n_elements(data.jd_ao)
 jac =  fltarr(2*nd+2*nd_ao, 13) ;13 elements of p
 pbinary =  p[5:11]
 rhotheta = binary_position(p[5:11], data.jd,  deriv = deriv)
 ;Now re-organise the Jacobian (see below for non-derivatives)
 dRAdrho = rebin(sin(rhotheta[*, 1]*!pi/180.0),nd,7)
 dRAdth  = rebin(rhotheta[*, 0]*!pi/180.*cos(rhotheta[*, 1]*!pi/180.0), nd,  7)
 jac[0:nd-1, 5:11] =dRAdrho*deriv[*, *, 0] + dRAdth*deriv[*, *, 1] 
 dDecdth = rebin(-rhotheta[*, 0]*!pi/180.*sin(rhotheta[*, 1]*!pi/180.0),nd,7)
 dDecdrho = rebin(cos(rhotheta[*, 1]*!pi/180.0),nd,7)
 jac[nd:2*nd-1, 5:11] =dDecdrho*deriv[*, *, 0] + dDecdth*deriv[*, *, 1]
 ;AO binary...
 pbinary[2] = p[12] ;Different semi-major axis for AO orbit.
 pbinary[4] = 180.0+pbinary[4] ;Secondary orbit, not primary...
 ;deriv_ao maps to p[5:11], except the second component which maps to p[12]
 rhotheta_ao =  binary_position(pbinary,  data.jd_ao, deriv = deriv_ao)
 ;Now for separation derivatives from AO
 jac[2*nd:2*nd+nd_ao-1, 5:6]  = deriv_ao[*, 0:1, 0]
 jac[2*nd:2*nd+nd_ao-1, 12]   = deriv_ao[*, 2, 0]
 jac[2*nd:2*nd+nd_ao-1, 8:11] = deriv_ao[*, 3:*, 0]
 jac[2*nd+nd_ao:2*nd+2*nd_ao-1, 5:6]  = deriv_ao[*, 0:1, 1]
 jac[2*nd+nd_ao:2*nd+2*nd_ao-1, 12]   = deriv_ao[*, 2, 1]
 jac[2*nd+nd_ao:2*nd+2*nd_ao-1, 8:11] = deriv_ao[*, 3:*, 1]
 radec_motion =  parallax_motion(data.jd, data.ra,  data.dec,  p[4])
 del_radec_motion = parallax_motion(data.jd, data.ra,  data.dec,  p[4]+0.01) - radec_motion
 jac[0:nd-1, 4] =  1e2*del_radec_motion[*, 0]
 jac[nd:2*nd-1, 4] =  1e2*del_radec_motion[*, 1]
 ra_plx_pm =  reform(radec_motion[*, 0]) + p[0] + p[2]*(data.jd-data.jd[0])/yr
 dec_plx_pm =  reform(radec_motion[*, 1]) + p[1] + p[3]*(data.jd-data.jd[0])/yr
 ;Some outputs to see the binary motion only
 bin_ra =  data.raoff - ra_plx_pm
 bin_dec =  data.decoff - dec_plx_pm
 modraoff =  ra_plx_pm + rhotheta[*, 0]*sin(rhotheta[*, 1]*!pi/180.0)
 moddecoff = dec_plx_pm + rhotheta[*, 0]*cos(rhotheta[*, 1]*!pi/180.0)
 jac[0:nd-1, 0] = 1.
 jac[nd:2*nd-1, 1] = 1.
 jac[0:nd-1, 2] = (data.jd-data.jd[0])/yr
 jac[nd:2*nd-1, 3] = (data.jd-data.jd[0])/yr
 retvect =  [(data.raoff - modraoff)/data.ra_err, (data.decoff - moddecoff)/data.dec_err, $
            (data.sep_ao-rhotheta_ao[*, 0])/data.sep_err, mod360(data.pa_ao-rhotheta_ao[*, 1])/data.pa_err]
 ;Last, the Jacobian has to be normalized for the data errors.
 for i = 0, 12 do jac[*, i] /= [data.ra_err, data.dec_err, data.sep_err, data.pa_err]
 return,  retvect
end
