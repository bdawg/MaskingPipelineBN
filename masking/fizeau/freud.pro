;######## RELEASE CLIP HERE ######## 
;+
; freud (head)
;
; freud will read the header passed to it, returning a structure with
; useful information extracted

; Input:  header     - the fits header file to be searched
; Returns: output     - a structure with everything you wanted to know
; 
; created                                                    PGT 01Dec03
; hacked-up to be simpler for LWS                            MJI 11Jan06

function freud,head, olog = olog,  ix = ix

; *** not tested yet ***

s='' &  i=0 &  f=float(0.0) &  d=double(0.0)

headinfo={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,optic_cfg:s, $
          lyot:s,grism:s,source_name:s,utc:s,date:s,jd:d,elevation:f,airmass:f,pa:f, $
          ra:d,dec:d,equinox:f, mask:s,  raoff:f, decoff:f,  del_elev:f,  del_pa:f}
;Allow a dummy headinfo to be returned.
if (head[0] eq '') then return,  headinfo

; Firstly figure out what sort of data we have ... 
iname=sxpar(head,'INSTRUME',count=count)
if(count eq 0) then iname=sxpar(head,'CURRINST',count=count)  ; this is for NIRC2

; NOTE need to add ra, dec, equinox to all but nirc2

if(strpos(iname,'NIRC') ne -1 and strpos(iname,'NIRC2') eq -1) then begin
  headinfo.instrument='NIRC'
  headinfo.nax1=sxpar(head,'NAXIS1')
  headinfo.nax2=sxpar(head,'NAXIS2')
  headinfo.nax3=sxpar(head,'NAXIS3')
  headinfo.t_int=sxpar(head,'TINT')
  headinfo.coadd=sxpar(head,'M56COAD0')
  headinfo.filter=strtrim(sxpar(head,'FILTER'),2)
  headinfo.slit=strtrim(sxpar(head,'SLTNAME'),2)
  headinfo.source_name=strtrim(sxpar(head,'TARGNAME'),2)
  headinfo.utc=strtrim(sxpar(head,'UTC'),2)
  headinfo.date=strtrim(sxpar(head,'DATE-OBS'),2)
  headinfo.jd=sxpar(head,'MJD-OBS')
  headinfo.elevation=sxpar(head,'EL')
  headinfo.airmass=sxpar(head,'AIRMASS')
  if(headinfo.jd lt 50500.) then $       ; This selects Jan97 & before where
     headinfo.pa=sxpar(head,'PARANG') $  ; U242N/U379N had different header keyword
     else headinfo.pa=sxpar(head,'PARANTEL')
endif

if(strpos(iname,'LWS') ne -1) then begin
  ;Straight to olog...
  olog.instrument[ix] = 'LWS'
  olog.nax1[ix]=sxpar(head,'NAXIS1')
  olog.nax2[ix]=sxpar(head,'NAXIS2')
  olog.nax3[ix]=sxpar(head,'NAXIS3')
  olog.t_int[ix]=sxpar(head,'FRMTIME')*sxpar(head,'FRMCOADD')
  olog.coadd[ix]=sxpar(head,'FRMCOADD') ;add 'CHPCOADD' ?
  olog.filter[ix]=strtrim(sxpar(head,'FILNAME'),2)
  olog.optic_cfg[ix]=strtrim(sxpar(head,'GRANAME'),2)
  olog.source_name[ix]=strtrim(sxpar(head,'TARGNAME'),2)
  olog.utc[ix]=strtrim(sxpar(head,'UTC'),2)
  olog.date[ix]=strtrim(sxpar(head,'DATE-OBS'),2)
  olog.jd[ix]=sxpar(head,'MJD-OBS')
  olog.elevation[ix]=sxpar(head,'EL')
  olog.airmass[ix]=sxpar(head,'AIRMASS')
  olog.pa[ix]= (sxpar(head,'PARANTEL') + sxpar(head,'ROTPPOSN') - 69.0) mod 360.0
  headinfo.instrument = 'LWS'
endif

if(strpos(iname,'TReCS') ne -1) then begin
  ;Straight to olog...
  ;!!! We need WINDOW for the pixel scale, and both filter wheels...
  headinfo.instrument = 'TRECS'
  olog.instrument[ix] = 'TRECS'
  olog.nax1[ix]=sxpar(head,'NAXIS1')
  olog.nax2[ix]=sxpar(head,'NAXIS2')
  olog.nax3[ix]=sxpar(head,'NAXIS3')
  olog.source_name[ix]=strtrim(sxpar(head,'OBJECT'),2)
  olog.utc[ix]=strtrim(sxpar(head,'TIME-OBS'),2)
  olog.date[ix]=strtrim(sxpar(head,'DATE-OBS'),2)
  olog.filter[ix]=strtrim(sxpar(head,'FILTER1'),2)
  olog.elevation[ix]=sxpar(head,'ELEVATIO')
  olog.airmass[ix]=sxpar(head,'AIRMASS')
  olog.t_int[ix]=sxpar(head,'OBJTIME')
  olog.coadd[ix]=sxpar(head,'FRMCOADD') ;add 'CHPCOADD' ?
  olog.jd[ix]=sxpar(head,'MJD-OBS')
  az =  sxpar(head,'AZIMUTH')
 ; lat =  -(30.+14./60.+16.8/3600.)
 ; altaz2hadec_rot,  olog.elevation[ix],  az,  lat,  ha,  dec,  par
 ; olog.pa[ix]= sxpar(head,'PA') + par ;!!! sign etc not done here !!!
  getrot,head,rot,cdelt
  print,  cdelt
  olog.pa[ix]= -rot
  maskname =  sxpar(head,'LYOT')
  case maskname of  
    'Ciardi  ' :  olog.mask =  't7'
     else     :  begin
     	    	    print,"NO MASK IN DATA"
		    stop
    	         endcase
  endcase
  olog.optic_cfg[ix]=strtrim(sxpar(head,'FILTER2'),2) ;!!! Check this...
endif

if(strpos(iname,'PHARO') ne -1) then begin
stop
  headinfo.instrument  = 'PHARO'
  olog.instrument[ix]  = 'PHARO'
  olog.ra[ix]          = strtrim(sxpar(head,'CRVAL1'),2)
  olog.dec[ix]         = strtrim(sxpar(head,'CRVAL2'),2)
  olog.ha[ix]          = strtrim(sxpar(head,  'HOURANGL'), 2)
  olog.filter[ix]      = strtrim(sxpar(head, 'FILTER'),2)
  olog.slit[ix]        = strtrim(sxpar(head, 'SLIT'),2)
  olog.optic_cfg[ix]   = strtrim(sxpar(head, 'CAROUSEL'),2)
  olog.lyot[ix]        = strtrim(sxpar(head, 'LYOT'),2)
  olog.grism[ix]       = strtrim(sxpar(head, 'GRISM'),2)
  olog.source_name[ix] = strtrim(sxpar(head, 'OBJECT'),2)
  olog.utc[ix]         = strtrim(sxpar(head, 'TIME-OBS'),2)
  olog.date[ix]        = strtrim(sxpar(head, 'DATE-OBS'),2)
  olog.nax1[ix]        = sxpar(head, 'NAXIS1')
  olog.nax2[ix]        = sxpar(head, 'NAXIS2')
  olog.nax3[ix]        = sxpar(head, 'NAXIS3')
  olog.equinox[ix]     = sxpar(head, 'EQUINOX')
  olog.airmass[ix]     = sxpar(head, 'AIR_MASS')
  olog.t_int[ix]       = sxpar(head, 'T_INT') * 1e-3
  ;From Stanmire Metchev's document:
  ; www.astro.caltech.edu/palomar/200inch/palao/Pharo/pharo_plate_scale.pdf
  ; and Mike's checks with Xi Cep from the first run.
;  olog.pa[ix]          = sxpar(head, 'CR_ANGLE') + 25.45
  olog.pa[ix]          = sxpar(head, 'CR_ANGLE') + 25.08 ;New Metchev value
  ; #### reads,(strmid(olog.date[0],0,4),thisyear)  ; Is it after P3K install?
  ; #### if(thisyear ge 2012.) then begin 
  ; ####  olog.pa[ix]          = (sxpar(head, 'CR_ANGLE')  -140.) mod 360 ;New P3K value BUT NOTE MUST ALSO FLIP X !!!
  ; #### endif 
  olog.coadd[ix]       = 1
  olog.raoff[ix]       =  sxpar(head, 'RA_OFFS')
  olog.decoff[ix]       =  sxpar(head, 'DEC_OFFS')
;  headinfo.raoff = sxpar(head, 'RAOFF')
;  headinfo.decoff = sxpar(head, 'DECOFF')

  maskname = sxpar(head, 'LYOT')
  case maskname of
      '18H     ' : olog.mask = 'p18'
      '9H - S2 ' : olog.mask = 'p9'
      '9H      ' : olog.mask = 'p9'
      else       : olog.mask = ''
  endcase
endif

if(strpos(iname,'NIRC2') ne -1) then begin
  headinfo.instrument='NIRC2'
  headinfo.nax1=sxpar(head,'NAXIS1')
  headinfo.nax2=sxpar(head,'NAXIS2')
  ; headinfo.nax3=sxpar(head,'NAXIS3')
  headinfo.t_int=sxpar(head,'ITIME')
  headinfo.coadd=sxpar(head,'COADDS')
  headinfo.filter=strtrim(sxpar(head,'FILTER'),2)
  headinfo.slit=strtrim(sxpar(head,'SLITNAME'),2)
  headinfo.lyot=strtrim(sxpar(head,'PMSNAME'),2)
  headinfo.optic_cfg=strtrim(sxpar(head,'CAMNAME'),2)
  headinfo.source_name=strtrim(sxpar(head,'TARGNAME'),2)
  headinfo.utc=strtrim(sxpar(head,'UTC'),2)
  headinfo.date=strtrim(sxpar(head,'DATE-OBS'),2)
  headinfo.jd=sxpar(head,'MJD-OBS')
  headinfo.elevation=sxpar(head,'EL')
  headinfo.airmass=sxpar(head,'AIRMASS')
  headinfo.pa=360.+sxpar(head,'PARANG')+sxpar(head,'ROTPPOSN') $
                  -sxpar(head,'EL')-sxpar(head,'INSTANGL')   
                                      ; THIS formula not checked yet!
  headinfo.ra=sxpar(head,'RA')
  headinfo.dec=sxpar(head,'DEC')
  headinfo.equinox=sxpar(head,'EQUINOX')
  headinfo.raoff = sxpar(head, 'RAOFF')
  headinfo.decoff = sxpar(head, 'DECOFF')
endif

; a bit hacked for now. Need more intelligence to work out filters (won't do for J)
if(strpos(iname,'CONICA') ne -1) then begin
  headinfo.instrument='CONICA'
  headinfo.nax1=sxpar_conica(head,'NAXIS1')
  headinfo.nax2=sxpar_conica(head,'NAXIS2')
  nax=sxpar_conica(head,'NAXIS')
  if(nax gt 2) then headinfo.nax3=sxpar_conica(head,'NAXIS3') 
  headinfo.t_int=sxpar_conica(head,'SODETDIT') ; HIERARCH ESO DET DIT
  ;headinfo.coadd=sxpar_conica(head,'COADDS')
  ; "filter" is tricky. Usually one of two wheels, but I have not yet covered case for J!
  wheel5=strtrim(sxpar_conica(head,'SOPTI5ID'),2) ; HIERARCH ESO INS OPTI6 ID
  wheel6=strtrim(sxpar_conica(head,'SOPTI6ID'),2) ; HIERARCH ESO INS OPTI6 ID
  if( (strpos(wheel5,"empty"))[0] eq -1) then headinfo.filter=wheel5 else headinfo.filter=wheel6
  ;headinfo.slit=strtrim(sxpar_conica(head,'SLITNAME'),2)
  headinfo.mask=strtrim(sxpar_conica(head,'SOPTI3ID'),2) ; HIERARCH ESO INS OPTI3 ID
  headinfo.lyot=strtrim(sxpar_conica(head,'SOPTI7ID'),2) ; ok, *NOT* lyot - this is the camera S13 or L27
  headinfo.optic_cfg=strtrim(sxpar_conica(head,'SOPTI4ID'),2) ; HIERARCH ESO INS OPTI4 ID - eg 'Wollaston_45'
  headinfo.source_name=strtrim(sxpar_conica(head,'OBJECT'),2) 
  if( (strpos(headinfo.source_name,"name not set"))[0] ne -1) then $
          headinfo.source_name=strtrim(sxpar_conica(head,'OOBSNAME'),2)
  headinfo.utc=strtrim(sxpar_conica(head,'UTC'),2)
  headinfo.date=strtrim(sxpar_conica(head,'DATE-OBS'),2)
  headinfo.jd=sxpar_conica(head,'MJD-OBS')
  headinfo.elevation=sxpar_conica(head,'EL')
  headinfo.airmass=sxpar_conica(head,'AIRMASS')
  ; Now work out sky PA of data cube. Not completely straightforward... why no start/end for Alt??
  rotstart=sxpar_conica(head,'ROTSTART')    ; HIERARCH ESO ADA ABSROT START
  rotend  =sxpar_conica(head,'BSROTEND')    ; HIERARCH ESO ADA ABSROT END
  pastart =sxpar_conica(head,'ANGSTART')    ; HIERARCH ESO TEL PARANG START
  paend   =sxpar_conica(head,'ARANGEND')    ; HIERARCH ESO TEL PARANG END
  alt     =sxpar_conica(head,'SOTELALT')    ; HIERARCH ESO TEL ALT
  instrument_offset= -0.55
  headinfo.pa=(rotstart+rotend)/2.+alt-(180.-(pastart+paend)/2.) + instrument_offset
  headinfo.ra=sxpar_conica(head,'RA')
  headinfo.dec=sxpar_conica(head,'DEC')
  headinfo.equinox=sxpar_conica(head,'GEQUINOX') ; HIERARCH ESO TEL TARG EQUINOX
  ;headinfo.raoff = sxpar_conica(head, 'RAOFF')
  ;headinfo.decoff = sxpar_conica(head, 'DECOFF')
endif

if(headinfo.instrument eq '') then $
  print,'*** ERROR - Freud could not understand the Header'

return,headinfo

end
 

