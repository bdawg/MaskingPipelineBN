;This function takes a wavelength vector (in microns) and a flux
;vector in Jy and returns a vector of magnitudes. It linearly
;interpolates between the bands, and uses values from Allen's
;Astrophysicsl quantities
;The option of filter = 'J'|'H'|'K'|'V'|'R'|'I' just finds the mean
;flux in the filter bandwidth, as defined roughly by Scholz model files
;and Allen.

function jy2mag, wavelengths, flux, filter=filter

;Values for a 0 mag star. Note that I've doctored the B filter a
;little for cubic spline interpolation, maitaining the mean flux
;within the filter.
filt_l =  [0.36,0.41,0.47,0.55, 0.7,0.9, 1.22,1.63,2.19,3.45,  5,  10.2,21.0]
filt_jy = [1840,4240,4320,3832,2842,2246,1570,1020,636,313,182.5,42.7,10.0];J was 1724
zmag_flux = interpol(filt_jy,filt_l,wavelengths, /spline)

if keyword_set(filter) then begin
 w = [-1]
 if (filter eq 'J') then w = where(wavelengths ge 1.1  and wavelengths le 1.4 )
 if (filter eq 'H') then w = where(wavelengths ge 1.483 and wavelengths le 1.783)
 if (filter eq 'K') then w = where(wavelengths ge 1.995 and wavelengths le 2.395)
 if (filter eq 'L') then w = where(wavelengths ge 3.04 and wavelengths le 3.84)
 if (filter eq 'L*') then w = where(wavelengths ge 3.5 and wavelengths le 4.1)
 if (filter eq 'V') then w = where(wavelengths ge 0.505 and wavelengths le 0.595)
 if (filter eq 'R') then w = where(wavelengths ge 0.595 and wavelengths le 0.805)
 if (filter eq 'I') then w = where(wavelengths ge 0.8 and wavelengths le 1.0)
 if (filter eq 'B') then w = where(wavelengths ge 0.39 and wavelengths le 0.49)
 if (filter eq '78') then w = where(wavelengths ge 0.7773 and wavelengths le 0.7863)
 if (filter eq '87') then w = where(wavelengths ge 0.8736 and wavelengths le 0.8818)
 if (filter eq '70') then w = where(wavelengths ge 0.695 and wavelengths le 0.705)
 if (filter eq '71') then w = where(wavelengths ge 0.710 and wavelengths le 0.720)
 if (filter eq '905 box') then w = where(wavelengths ge 0.878  and wavelengths le 0.930)
 if (filter eq '1040 narrow') then w = where(wavelengths ge 1.035  and wavelengths le 1.045)
 if (w[0] eq -1) then begin
  print, 'Error: No input wavelengths within filter...'
  return, 0
 endif
 nw = n_elements(w)
 mns = 0.5*(flux[w[0:nw-2]] + flux[w[1:nw-1]])
 diffs = wavelengths[w[1:nw-1]] - wavelengths[w[0:nw-2]]
 mnl = 0.5*(wavelengths[w[1:nw-1]] +wavelengths[w[0:nw-2]])
 weights = replicate(1.0,n_elements(mns)) ;an array of the same size...
 if (filter eq 'L') then weights = interpol([0,0.5,0.85,0.86,0.19,0],[3.04,3.24,3.48,3.68,3.76,3.84],mnl)
 if (filter eq 'V') then weights = interpol([0,0.78,1.0,0.898,0.359,0.08,0],[0.48,0.51,0.53,0.55,0.6,0.64,0.68],mnl)
 mnflux = total(mns*weights*diffs/mnl^2)/total(diffs/mnl^2)
;Now do the same for zmag_flux
 zmag_flux = 0.5*(zmag_flux[w[0:nw-2]] + zmag_flux[w[1:nw-1]])
 mnzflux = total(zmag_flux*weights*diffs/mnl^2)/total(diffs/mnl^2)
 return, 2.5*alog10(mnzflux/mnflux)
endif

return, 2.5*alog10(zmag_flux/flux)

end
