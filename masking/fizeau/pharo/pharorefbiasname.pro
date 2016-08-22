;$Id: pharorefbiasname.pro,v 1.1.1.1 2005/12/19 05:15:05 mireland Exp $
;$log$


function pharorefbiasname, header
;
; take a header and make a standard refbias name
;
return, 'phrefbias.'+pharomode(header)+'.fits'

end

;prefix = sxpar(header, 'OBJECT')
;return, 'pharo.'+prefix+pharomode(header)+'.fits'
