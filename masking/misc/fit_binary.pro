;This quick procedure is used by mpfit2dfunc to fit a binary
;separation and position angle.
;b:     array of baseline lengths in wavelengths
;theta: array of baseline pas
;parms: 4 parameters, separation, position angle, brightness ratio,
;  system correlation...OR add diam or primary and diam of secondary
;p[0]: separation (mas)
;p[1]: PA (degrees)
;p[2]: brightness ratio (star 1/star 2)
;p[3]: system correlation
;p[4]: diam 1
;p[5]: diam 2

function fit_binary, b, theta, parm_in
  
p = parm_in
if (n_elements(p) eq 4) then p = [p,mas2rad(0.001),mas2rad(0.001)]
 ro = mas2rad(p[0])
 alpha = p[1]*!pi/180.0
; minv2 = ( (1-p[2])/(1+p[2]) )^2 * p[3]
; maxv2 = p[3]
 uniform_disk, b, [mas2rad(p[4]),1.0], v1
 uniform_disk, b, [mas2rad(p[5]),1.0], v2
 result = p[2]^2*v1^2 + v2^2  + 2*v1*v2*p[2]*cos(2*!pi*b*ro*(cos(theta)*cos(alpha)+sin(theta)*sin(alpha)))
 
 return, p[3]*result/(1+p[2])^2

end
