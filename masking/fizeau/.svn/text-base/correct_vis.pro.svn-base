;######## RELEASE CLIP HERE ######## 
;This function corrects visibility (or V^2) for seeing and
;windshake. Use this on source and cal before division.
;Returns a multiplicative correction to the 
;Input parameters:
; vis:   vis or V^2
; avar:  variance in a single frame of data due to atmosphere
; corr_const: Correction constant. 0.4 is good for V for NIRC
; experiment.

function correct_vis, vis, avar, u,v, corr_const, err_avar=err_avar, err_vis=err_vis, nows=nows

w = where((vis eq 0) or (avar eq 0))
if (w[0] ne -1) then begin
 print, 'Cannot continue - divide by zero imminent...'
 stop
endif
nf = 100.0 ;Not really important...
normvar = avar/vis^2
w0 = where(u^2+v^2 gt 0)
parinfo = replicate({fixed:0}, 5)
if (keyword_set(nows)) then parinfo[3:4].fixed = 1
if (keyword_set(err_avar) eq 0) then err_avar =  normvar*2.0/sqrt(nf)
if (keyword_set(err_vis)) then err_avar = abs(normvar)*sqrt(err_avar^2/avar^2 + 2.0*err_vis^2/vis^2)
try = mpfit2dfun('v2varfunc', u[w0], v[w0], normvar[w0] , err_avar[w0], $
                 [0,0,0,0,0], /quiet, yfit=yfit, perror=perror)
correction = vis ;Make correction an array identical to vis, filled with 1.0.
correction[*] = 1.0
correction[w0] = exp(yfit*corr_const)
return, correction

end

