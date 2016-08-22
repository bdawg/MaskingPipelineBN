;
; Feb 19, 1997
; JDM
;
; angle_diff
;
;   Pass this routine two angles, or arrays of angles, and
;   this routine will return to you the Difference between them.
;
;  ** DEFAULT is in DEGREES! **
;
function angle_diff,angles1,angles2,radians=radians,help=help

if (keyword_set(help) eq 1) then begin
 print,' function angle_diff,angles1,angles2,radians=radians,help=help'
 print, '  DEFAULT is in DEGREES, else use /radians'
endif

a1=angles1/!radeg
a2=angles2/!radeg

if (keyword_set(radians) eq 1) then begin
  a1=a1*!radeg
  a2=a2*!radeg
endif

; a1,a2 in radians
;

;cosdiff=cos(a1)*cos(a2)+sin(a1)*sin(a2)
; sindiff=sin(a1)*cos(a2)-cos(a1)*sin(a2)
;diff=acos(cosdiff)
; diff=asin(sindiff)
diff=(mod360(!radeg*a1-!radeg*a2))
;diff=diff*!radeg
if (keyword_set(radians) eq 1) then begin
 diff=diff/!radeg
endif else $
; diff=mod360(diff)
return,diff
end



