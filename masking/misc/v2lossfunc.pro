;This predicts the V^2 loss from normalised v2 variance...

function v2lossfunc, seevar, wsvar, parms

 return, parms[0]*exp(-seevar/parms[1] - wsvar/parms[2])
   
end
