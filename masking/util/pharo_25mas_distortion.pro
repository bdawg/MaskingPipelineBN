pro pharo_25mas_distortion, x, y, ha, dec, cr_ang=cr_ang, crN=crN, $
	xdp=xdp, ydp=ydp, xp=xp, yp=yp, xc=xc, yc=yc, mag=mag, angle=angle

; pharo_25mas_distortion	- Produce the corrected (xdp, ydp) 
;			coordinates at a given location (x,y) on the 
;			PHARO 25mas/pix array, for given
;			hour angle (ha), declination (dec), and
;			Cassegrain ring angle (cr_ang).
;
; x, y		- Input pixel coordinates.
; ha		- Telescope hour angle (west is negative).
; dec		- Telescope declination.
; cr_ang	- Cassegrain ring angle (optional).
; crN		- Cassegrain ring angle at which celestial North and the
;		  detector y-axis are aligned (default=334.88degrees).
; xdp, ydp	- Output distortion-corrected coordinates.
; xp, yp	- Output coordinates before the beam tilt correction.
; xc, yc	- Origin of the (x,y) system; (default = 512,512).
; angle		- Orientation of the beam tilt (default=82.5 degrees).
; mag		- Beam tilt magnification along the angle+90 axis
;		  (default=1.0069).

if n_params(0) lt 4 then $
    message,'Usage: PHARO_25MAS_DISTORTION, x, y, hour_angle, dec'+$
	'[, cr_ang=cr_ang, crN=crN, xdp=xdp, ydp=ydp, xp=xp, yp=yp, xc=xc,'+$
	' mag=mag, angle=angle]'
if n_elements(xc) lt 1 then xc = 512.
if n_elements(yc) lt 1 then yc = 512.
if n_elements(mag) eq 0 then mag = 1.0069	; from May 3, 3005
if n_elements(angle) eq 0 then angle = 82.5	; from May 3, 2005
if n_elements(crN) eq 0 then crN= 334.88	; from May 3, (11?), 2005

; Convert to radians
pi = 3.141592653D
beamtilt_ang = angle/180.*pi
northPA = crN/180.*pi
if n_elements(cr_ang) gt 0 then cassring_ang = cr_ang/180.*pi $
else cassring_ang = northPA

x = x-xc
y = y-yc

; Set the xp coefficients in xp(x,y) = a0+a1*y+a2*y^2+a3*x+a5x^2
a = dblarr(6)
b = dblarr(3)
a1 = get_a1 (ha, dec)
a[0] = 0.01839547	; 300
a[1] = a1
a[2] = -4.510687e-06	; 300
a[3] = 0.9998296 	; 300
a[4] = 0
a[5] = -4.365845e-06	; 300

; Set the yp coefficients in yp(x,y) = b0+b1*y+b2*x
b[0] = 0.01932465	; 300
b[1] = 1.00001000	; 300
b2 = get_b2 (ha, dec)
b[2] = b2

; Correct for the detector distortion
xp = a[0] + a[1]*y + a[2]*y^2 + a[3]*x + a[4]*y*x + a[5]*x^2 
yp = b[0] + b[1]*y + b[2]*x

; Correct for the distortion due to the primary+secondary beams.
; Project along the (beamtilt_ang,beamtilt_ang+90) coordinate set and apply 
; the supplied magnification factor mag along the beamtilt_ang+90 axis.  
; Then de-project without de-magnifying.  The beam tilt angle relative
; to the detector also depends on the CR angle.  Tested - it seems that
; most of the distortion is unrelated to the primary+secondary, but comes
; from the positioning of the detector in the focal plane.
xr = xp*cos(beamtilt_ang)+yp*sin(beamtilt_ang)
yr = (-xp*sin(beamtilt_ang)+yp*cos(beamtilt_ang)) * mag
xdp = xr*cos(beamtilt_ang)-yr*sin(beamtilt_ang)
ydp = xr*sin(beamtilt_ang)+yr*cos(beamtilt_ang)

; Rotate so that N is exactly up
if (n_elements(cr_ang) gt 0) then begin
	theta = cassring_ang - northPA
	xdpN = xdp*cos(theta)-ydp*sin(theta)
	ydpN = xdp*sin(theta)+ydp*cos(theta)
	xdp = xdpN
	ydp = ydpN
endif

; The distortion-corrected coordinates
xdp += xc
ydp += yc

; Restore the input coordinates to their original values
x += xc
y += yc

end
