function mod360, angles

; this routine will makes an angle (in degrees) to be between
; -180 and 180 degrees.

;newangles=((((angles+180.d) mod 360.d)+360.d) mod 360.d) - 180.d
;; mod returns a number of the same sign, but we don't want
;; that.  The 10980 accomplishes this practically without needing two mods.
newangles = ( (angles + 10980.d) mod 360.d ) - 180.d

return,newangles
end



