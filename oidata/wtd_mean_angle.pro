; JDM Feb 26, 1997
; JDM 2005May30  using more standard method.

;
; This routine will Accept a bunch of ANGLES and optional Weights
; and will return the MEAN Angle and the derived STDEV of this MEAN.
;
; result( which is MEAN, SIGMA) = angle_stats(angles,weights=weights)
;  Result and angles in Degrees unless specified otherwise via /radians
;

function wtd_mean_angle,angles,err,radians=radians,werr 
; output wtd_mean with option WERR for error on mean

if (keyword_set(angles) eq 0) then begin
  print,'function angle_stats,angles,weights=weights
  print,' RETURNS  MEAN ANGLE, STDEV of MEAN.'
  print,' ANGLES IN DEGREES!!!!'
endif

if (keyword_set(radians) eq 1) then factor=1.0 else factor=!pi/180.

num=n_elements(angles)
if (keyword_set(err) eq 0) then err=replicate(1.0,num)

                rp_array=cos(factor*angles)
		ip_array=sin(factor*angles)
	  	rp_err=tan(factor*err)
                rp_mean=wtd_mean(rp_array,rp_err,rp_sigma2)
                ip_mean=wtd_mean(ip_array,rp_err,ip_sigma2)
                ri2at,rp_mean,ip_mean,amp,theta
                mean_theta=mod360(theta)
                mean_sigma=!radeg*atan( (rp_sigma2 > ip_sigma2 )/amp)
result=mean_theta*!pi/180. / factor
werr=mean_sigma  *  !pi/180. / factor
return,result
end




