;This program gives delta Dec, delta RA based on:
; JD, RA (degrees), Dec (degrees) and parallax
;jd = julday(0,0,2006,0,0,0) +findgen(365)
;ra = (20+43/60.+19.4/3600.)*15.0
;dec = 55+20/60. +52.0/3600.
;a = read_ascii('gj802jul06dat.txt')
;a =a.field1[*,0:25]
;raoffset =  a[1,*]-(a[0,*]-a[0,0])*2.37 
;decoffset = a[2,*]-(a[0,*]-a[0,0])*4.735
function parallax_motion,  jd,  ra,  dec,  parallax

jd =  reform(jd)
;Length of year from Wikipedia
yr = double(365.256363051)
;Earth's current tilt:
eps =  23.439281*!pi/180.0
;Parallax in radians
prad =  mas2rad(parallax)

;Ecliptic Longitude of the earth.
;Date of Vernal Equinox from:
;http://aa.usno.navy.mil/data/docs/EarthSeasons.html
;0 hours is at midnight in September.
;long_e =  jd - julday(3, 20, 2000, 7, 35, 00)
long_e =  jd - julday(9, 22, 2000, 17, 27, 00)
long_e =  ((long_e + 100.0*double(yr)) mod yr)/yr*2*!dpi

;Three equations from Wikipedia, to get the ecliptic lattitude and
;longitude of the star in radians
e1 =  cos(eps)*sin(dec*!dpi/180.0) - sin(ra*!dpi/180.0)*cos(dec*!dpi/180.0)*sin(eps)
e2 =  cos(ra*!dpi/180.0)*cos(dec*!dpi/180.0)
e3 =  sin(eps)*sin(dec*!dpi/180.0) + sin(ra*!dpi/180.0)*cos(dec*!dpi/180.0)*cos(eps)
long_s =  atan(e3, e2)
lat_s  =  atan(e1, e2/cos(long_s))

;Now find delta_long and delta_lat...
;NB The - sign on the next line is a mystery. Also, the simple
;sin and cos orbit should change to an elliptical one, by:
HELIO, jd, 3, HRAD, HLONG, HLAT
long_e =  hlong*!dpi/180.0 ;More accurate

delta_long =  -hrad*sin(long_e - long_s)*prad/cos(lat_s)
delta_lat  =  hrad*cos(long_e - long_s)*prad*sin(lat_s) 

;Temp code? Re-calculate old ra and dec using inverse Wikipedia equations
e1 =  sin(eps)*sin(long_s)*cos(lat_s) + cos(eps)*sin(lat_s)
e2 =  cos(long_s)*cos(lat_s)
e3 =  cos(eps)*sin(long_s)*cos(lat_s)-sin(eps)*sin(lat_s)
old_ra =  atan(e3, e2)*180.0/!dpi
if (old_ra[0] lt 0) then old_ra += 360.0
old_dec =  atan(e1, e2/cos(old_ra*!dpi/180.0))*180.0/!dpi

;Finally, use the three inverse Wikipedia equations to find the
;equatorial coordinates of the star, and the shifted equitorial coords
long_s +=  delta_long
lat_s  += delta_lat
e1 =  sin(eps)*sin(long_s)*cos(lat_s) + cos(eps)*sin(lat_s)
e2 =  cos(long_s)*cos(lat_s)
e3 =  cos(eps)*sin(long_s)*cos(lat_s)-sin(eps)*sin(lat_s)
new_ra =  atan(e3, e2)*180.0/!dpi
if (new_ra[0] lt 0) then new_ra += 360.0
new_dec =  atan(e1, e2/cos(new_ra*!dpi/180.0))*180.0/!dpi


return,  [[(new_ra-old_ra)*cos(old_dec[0]*!pi/180.)],  [new_dec-old_dec]]*3600000.0 ;In mas now.
;return,  [[new_ra-old_ra],  [new_dec-old_dec]]*3600000.0 ;In mas now.

end
