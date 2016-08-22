;This is a function to be fed to mpfit2dfun for fitting to normv2var
;to a 5-parameter windshake-seeing function

function v2varfunc, b_lengths, b_angles, parms

 return, parms[0] + parms[1]*b_lengths + parms[2]*b_lengths^2 + $
   parms[3]*b_lengths*abs(sin(b_angles)) + parms[4]*b_lengths^2*sin(b_angles)^2
   
end
