;This quick procedure is used by mpfit2dfunc to fit a binary
;separation and position angle.
;b:     array of baseline lengths
;theta: array of baseline pas
;parms: 4 parameters, separation, position angle, brightness ratio,
;  system correlation

function fit_binary2, b, theta, parms
  
 parms[2] = 0.58
 
 ro = mas2rad(parms[0])
 alpha = parms[1]*!pi/180.0
 minv2 = ( (1-parms[2])/(1+parms[2]) )^2 * parms[3]
 maxv2 = parms[3]
 return, (minv2+maxv2)*0.5 + 0.5*(maxv2-minv2)*cos(2*!pi*b*ro*(cos(theta)*cos(alpha)+sin(theta)*sin(alpha)))

end
