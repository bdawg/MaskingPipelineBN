; Does the job of Freud but for VAMPIRES data, when headers are stored
; in the extension 1 of the fits file (as opposed to Lacan, which used logfiles).
;
; Usage: hlog=jung(datadir,fitsfilename)

;; Email from Minowa-san:
;; I guess you would like to use ADI mode what we call, which fix the
;; rotation of the pupil
;; on the science camera detector.
;; If you use the ADI mode with PA=0 deg,
;; the horizontal axis (parallel to the AO188 optical bench) is the
;; telescope elevation direction.

;; As you wrote in you e-mail, if you want to keep the top of the
;; detector pointing to the vertical top direction,
;; PA of the ADI mode should be -90deg or 90deg (I cannot tell which sign
;; is correct, since it depends on the combination of the camera lens).

;; D_IMRPAP is the PA (or offset angle) of the ADI mode.
;; If you want to know the angle between the celestial north on the
;; camera image and the top of the detector,
;; you need to refer D_IMRPAD.  This angle is determined by a combination
;; of the parallactic angle and image rotator. Again, I cannot tell zero
;; point and sign of the angle as it depends on the camera.


function jung,head


s='' &  i=0 &  f=float(0.0) &  d=double(0.0)

headinfo={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,optic_cfg:s, $
          lyot:s,grism:s,source_name:s,utc:s,date:s,jd:d,elevation:f,airmass:f,pa:f, $
          ra:d,dec:d,equinox:f, mask:s,  raoff:f, decoff:f,  del_elev:f,  del_pa:f, $
          emgain:f, hwp:f, timingpattern:s, imgRotAng:f, imgRotPad:f, imgRotPap:f ,$
          ADCstagepos:f, ADCp1angle:f, ADCp2angle:f, azimuth:f, localtime:s, pqwp1:f, pqwp2:f} ;These last 2 lines are added for VAMPIRES

; What is the PA? It is IMRPAD plus some instrument offset. Store PAD
; as PA, and rely on adding instrument offset later in the process.
PA = sxpar(head,'PAD')


headinfo.instrument='VAMPIRES'
headinfo.filter=sxpar(head,'FILTER')
headinfo.utc=sxpar(head,'UTSTTIME')
headinfo.airmass=sxpar(head,'airmass')
headinfo.pa=PA
headinfo.ra=ten(sxpar(head,'RA'))
headinfo.dec=ten(sxpar(head,'DEC'))
headinfo.localtime=sxpar(head,'loctime')
headinfo.emgain=sxpar(head,'emgain')
headinfo.hwp=sxpar(head,'A0HWPVAL')
headinfo.mask=sxpar(head,'MASK')
headinfo.timingpattern=sxpar(head,'TIMING')
headinfo.imgRotAng=sxpar(head,'IMRA')
headinfo.imgRotPad=sxpar(head,'PAD')
headinfo.imgRotPap=sxpar(head,'PAP')
headinfo.pqwp1=sxpar(head,'QWP1')
headinfo.pqwp1=sxpar(head,'QWP2')

return,headinfo

end
