;;This exracts useful information from the headers.
dir = '~snert/conica/data/conica_Mar08/15Mar08/NACO_IMG_SCI076_'
suffix = '_DIT.fits.gz'
nf = 200
;-------
for i=1,nf do begin
 h = headfits(dir + string(i, format='(I04)')+suffix)
; print, sxpar_conica(h, 'INSGRPID')
 print, i, sxpar_conica(h, 'SOPTI4ID')
endfor
;HIERARCH ESO INS GRP NAME
end
