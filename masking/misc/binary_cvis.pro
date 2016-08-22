;This function returns the complex visibilities for a binary based on
;a set of parameters.
;u, v : u,v coordinates...
;parms: 3 parameters, separation, position angle, brightness ratio

function binary_cvis, u, v, parms
  
 ro = mas2rad(parms[0])
 alpha = parms[1]*!pi/180.0 ; converts to radians
 return, (1.0 + parms[2]*exp(complex(0,2*!pi*ro*(u*cos(alpha) + v*sin(alpha)))))/(1.0 + parms[2])

end
