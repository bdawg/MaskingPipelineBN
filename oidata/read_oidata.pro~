pro read_oidata, file, oiarray, oitarget,oiwavelength,$
   oivis,oivis2, oit3, inventory=inventory
;+
; NAME:
;     read_oidata
; PURPOSE:
;     To read FITS files saved in the OI-DATA format into IDL structures
;	for analysis
;
; CALLING SEQUENCE:
;     read_oidata, fits_file, oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
;                       or
;     read_oidata, fits_file, oiarray,oitarget,oiwavelength,dummy,oivis2
;
; INPUTS:
;     fits_file - Name of fits file to read. Must be in OI-DATA format. 
;		  See Modification History to see what Revision #s are 
;		  supported and for URL of the OI-DATA format definitions.
;
; OPTIONAL INPUTS:
;     /INVENTORY Print message displaying the inventory of binary tables found 
;  		 in the OI-FITS file
;
; OUTPUTS:
;
;      OIARRAY   Contains information on the interferometer array geometry.
; 		 IDL array of OIARRAY structures containing all the information
;		 in the OI_ARRAY Binary FITS Table if present in the OIDATA
;		 FITS file.  There is one entry for each row in the binary 
;		 table (one row for each telescope) and the common header 
;		 keywords are repeated in each row.  Different OIARRAY tables
;                are distinguished using the ARRNAME (and EXTVER) keyword.
;      OITARGET	 Contains information on the observed targets. 
;                IDL array of OITARGET structures containing all the information
;                in the OI_TARGET Binary FITS Table if present in the OIDATA
;                FITS file.  There is one entry for each row in the binary 
;                table (one row for each target) and the common header 
;                keywords are repeated in each row. Only one OITARGET table
;                allowed.
;   OIWAVELENGTH Contains information on the observing bandpasses.
;                IDL array of OIWAVELENGTH structures containing all the 
;                information in the OI_WAVELENGTH Binary FITS Table if present
;	         in the OIDATA FITS file.  There is one entry for each instance 
;		 of the OIWAVELENGTH table and we note that EFF_WAVE/EFF_BAND 
;		 are pointers, not arrays themselves. Different OIWAVELENGTH 
;                tables are differentiated also by different INSNAME keywords.
;      OIVIS	 Contains complex visibility measurements.
;		 IDL array of OIVIS structures containing all the 
;                information in the OI_VIS Binary FITS Table if present
;                in the OIDATA FITS file.  There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row. 
;                Note that the structures contain pointers to the data, not
;                the data themselves. Also, different instances of each OIVIS
;                table are differentiated using combination of ARRNAME/INSNAME
;                and DATE-OBS.
;      OIVIS2    Contains visibility-squared measurements.
;                IDL array of OIVIS2 structures containing all the 
;                information in the OI_VIS2 Binary FITS Table if present
;                in the OIDATA FITS file.  There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row.
;                Note that the structures contain pointers to the data, not
;                the data themselves. Also, different instances of each OIVIS
;                table are differentiated using combination of ARRNAME/INSNAME
;                and DATE-OBS.
;      OIT3      Contains measurements of the bispectrum (triple amp and phase).
;                IDL array of OIT3 structures containing all the 
;                information in the OI_T3 Binary FITS Table if present
;                in the OIDATA FITS file.  There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row.
;                Note that the structures contain pointers to the data, not
;                the data themselves. Also, different instances of each OIVIS
;                table are differentiated using combination of ARRNAME/INSNAME
;                and DATE-OBS.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      There are currently 6 different binary tables defined in the OIDATA
;      format.  All tables are not required to be present to be a valid OIDATA 
;      FITS file.  READ_OIDATA will read the file and check for the
;      presence of the 6 allowed tables.  The data from each table is then read
;      into IDL structure arrays.  The structures are defined in procedures
;      DEFINE_OI*.PRO elsewhere in the IOI Library. 
;
;      The Data Exchange Standard for Optical (Visible/IR) Interferometry
;      is being maintained by the IAU Working group on Optical Interferometry
;      and the current version of the standard can be found at :
;            http://www.mrao.cam.ac.uk/~jsy1001/exchange/
;
; EXAMPLE:
;       Read in the sample OI-DATA FITS file distributed by IAU-WG:
;
;       IDL> READ_OIDATA, 'testdata.fits', oiarray,oitarget,oiwavelength,$
;		oivis, oivis2,oit3, /inventory
;	This file Satisfies the requirements of the OI_DATA format
;	Inventory:
;	  OI_ARRAY:             1
;	  OI_TARGET:            1
;	  OI_WAVELENGTH:        1
;	  OI_VIS:               1
;	  OI_VIS2:              1
;	  OI_T3:                1
;	 Unknown Tables:        0
;	IDL> print, oitarget.target
;		==>  alp_aur
;	IDL> print,oivis2(0).time,*oivis2(0).vis2data,*oivis2(0).vis2err
; 		==>        82810.000      0.67700000     0.064000000
;       IDL> print,*oivis(0).flag,':',string(*oivis(0).flag)
;		==> 	   70:F
;
;
; PROCEDURES USED:
;	This routine is part of the IDL Optical Interferometry (IOI) Library
;	 (more information at www.astro.lsa.umich.edu/~monnier).
;       The routines of the IOI Library generically 
;       require an up-to-date Astrolib library in order to read and write binary
;       FITS tables, among other things. The IDL Astronomy User's Library can
;       be found at http://idlastro.gsfc.nasa.gov/homepage.html.
;	
;
; MODIFICATION HISTORY:
;     v0.0 2002Jul04    J.D. Monnier    Written for OI-DATA Revision 0
;     v0.1 2003Feb07    J.D. monnier    Updated to reflect 2002Nov26 update
;          of the draft standard (Still officially Revision 0)
;          Adopted use of pointers to variable-length arrays (NWAVE) and
;          included support of multiple instances of OI-Tables.
;     v0.2 2003Feb13  	J.D. monnier	Updated to reflect oi-TARGET changes
;     v0.3 2003Feb18    J.D. Monnier	Use concat_oitable() to make routines
;	   compatibile back to IDL 5.1 at least (was only good for >5.4).
;     v0.4 2003Jul03    JDM small change to destroy oi* variable
;		when called so not to mix previous values with new ones
;
;     To do: A.  Add Revision Checking (once multiple revisions exist)  
;	     B.  Error checking if all columns headings and keywords are not used
;                (should fill in default dummy values)
;
;  *** Please write monnier@umich.edu with bug reports. ***
;-

; First destroy output variables
delvarx, oiarray, oiwavelength, oitarget,oivis,oivis2,oit3


; READING OI DATA
;--------------------
; Determine which binary tables are present and load units.
; units = [oiarray,oitarget,oiwavelength,oivis,oivis2,oit3]
; if unit lt 0 then the table does not exist in OIDATA fits file.
; According to IAU definition Revision 0, each file must contain a single
; OI_TARGET table, at least one OI_WAVELENGTH table, and AT LEAST ONE 
; data table (OI_VIS, OI_VIS2, or OI_T3).  In addition, the format
; supports multiple instances of OI_ARRAY, OI_WAVELENGTH, OI_VIS, OI_VIS2, and
; OI_T3 tables, and this reader purports to support this.

; Read HDU (Header and Data Unit)
errmsg=''
fxread,file,data,head_hdu,errmsg=errmsg
if (strmid(errmsg,0,5) eq 'Error') then begin
 print, errmsg + "   < RETURNING > "
 return
endif


; Find the available binary extensions
;
; Loop through all the extensions.
exten=1 ;total extensions
extver_array=0
extver_target =0
extver_wavelength=0
extver_vis=0
extver_vis2=0
extver_t3=0
unknown_table=0

errmsg=''
fxbopen, unit,file,exten,header,errmsg=errmsg

while (errmsg eq '')  do begin 

extname = strtrim(fxpar(header,'EXTNAME'),2 )
;print,'extname: ',extname
;help,extname
if (extname eq 'OI_ARRAY') then begin   ; Read in OI_ARRAY
 extver_array=extver_array+1
 head_oiarray=header
 ntel=fxpar(head_oiarray,'NAXIS2')
 oi_revn =  fxpar(head_oiarray,'OI_REVN')
 arrname =  fxpar(head_oiarray,'ARRNAME')
 frame =  fxpar(head_oiarray,'FRAME')
 arrayx =  fxpar(head_oiarray,'ARRAYX')
 arrayy =  fxpar(head_oiarray,'ARRAYY')
 arrayz =  fxpar(head_oiarray,'ARRAYZ')

 fxbread,unit,tel_name,'TEL_NAME'
 fxbread,unit,sta_name,'STA_NAME'
 fxbread,unit,sta_index,'STA_INDEX'
 fxbread,unit,diameter,'DIAMETER'
 fxbread,unit,staxyz,'STAXYZ'

; Define and Fill oi_array structure
  define_oiarray,oiarray
  oiarray=replicate(oiarray,ntel)
  oiarray(*).extver = extver_array
  oiarray(*).oi_revn=oi_revn
  oiarray(*).arrname =arrname
  oiarray(*).frame = frame
  oiarray(*).arrayx=arrayx
  oiarray(*).arrayy=arrayy
  oiarray(*).arrayz=arrayz
  oiarray.tel_name=tel_name
  oiarray.sta_name=sta_name
  oiarray.sta_index=sta_index
  oiarray.diameter=diameter
  oiarray.staxyz=staxyz

; If previous oiarray exists, then append; otherwise create it
  if (n_elements(oiarray_full) eq 0) then oiarray_full=oiarray $
   else oiarray_full=concat_oitable(oiarray_full,oiarray)

endif ; OI_ARRAY

if (extname eq 'OI_TARGET') then begin   ; Read in OI_WAVELENGTH
 extver_target=extver_target+1
 head_oitarget=header
 
 oi_revn =  fxpar(head_oitarget,'OI_REVN')
 ntargets=fxpar(head_oitarget,'NAXIS2')
 fxbread,unit,target_id,'TARGET_ID'
 fxbread,unit,target,'TARGET'
 fxbread,unit,raep0,'RAEP0'
 fxbread,unit,decep0,'DECEP0'
 fxbread,unit,equinox,'EQUINOX'
 fxbread,unit,RA_ERR,'RA_ERR'
 fxbread,unit,DEC_ERR,'DEC_ERR'
 fxbread,unit,SYSVEL,'SYSVEL'
 fxbread,unit,VELTYP,'VELTYP'
 fxbread,unit,VELDEF,'VELDEF'
 fxbread,unit,PMRA,'PMRA'
 fxbread,unit,PMDEC,'PMDEC'
 fxbread,unit,PMRA_ERR,'PMRA_ERR'
 fxbread,unit,PMDEC_ERR,'PMDEC_ERR'
 fxbread,unit,PARALLAX,'PARALLAX'
 fxbread,unit,PARA_ERR,'PARA_ERR'
 fxbread,unit,SPECTYP,'SPECTYP'

; Fill oi_array structure
 define_oitarget,oitarget
 oitarget=replicate(oitarget,ntargets)
  oitarget(*).oi_revn=oi_revn
  oitarget.target_id =target_id
  oitarget.target = target
  oitarget.raep0=raep0
  oitarget.decep0=decep0
  oitarget.equinox=equinox
  oitarget.ra_err=ra_err
  oitarget.dec_err=dec_err
  oitarget.sysvel=sysvel
  oitarget.veltyp=veltyp
  oitarget.veldef=veldef
  oitarget.pmra=pmra
  oitarget.pmdec=pmdec
  oitarget.pmra_err=pmra_err
  oitarget.pmdec_err=pmdec_err
  oitarget.parallax=parallax
  oitarget.para_err=para_err
  oitarget.spectyp=spectyp


; If previous oiarray exists, then append; otherwise create it
  if (n_elements(oitarget_full) eq 0) then oitarget_full=oitarget $
   else oitarget_full=concat_oitable(oitarget_full,oitarget)

endif ; OI_TARGET


if (extname eq 'OI_WAVELENGTH') then begin   ; Read in OI_WAVELENGTH
 extver_wavelength=extver_wavelength+1
 head_oiwavelength=header

 nwave=fxpar(head_oiwavelength,'NAXIS2')
 oi_revn =  fxpar(head_oiwavelength,'OI_REVN')
 insname =  fxpar(head_oiwavelength,'INSNAME')
 fxbread,unit,eff_wave,'EFF_WAVE'
 fxbread,unit,eff_band,'EFF_BAND'

; Current scheme: the array of wavelengths are stored in a vector in the
; the oi_wavelength structure, not as rows (which is the OI_WAVELENGTH format standard) 

 define_oiwavelength,oiwavelength,nwave=nwave
 oiwavelength.extver = extver_wavelength
 oiwavelength.nwave  = nwave ; implicitly also done in define_call above
 oiwavelength.oi_revn = oi_revn
 oiwavelength.insname = insname
 oiwavelength.eff_wave=ptr_new(eff_wave)
 oiwavelength.eff_band=ptr_new(eff_band)


; If previous array exists, then append; otherwise create it
  if (n_elements(oiwavelength_full) eq 0) then oiwavelength_full=oiwavelength $
   else oiwavelength_full=concat_oitable(oiwavelength_full,oiwavelength)

endif ; OI_WAVELENGTH


if (extname eq 'OI_VIS') then begin   ; Read in OI_VIS
 extver_vis=extver_vis+1
 head_oivis=header

 nvis=fxpar(head_oivis,'NAXIS2')
 oi_revn =  fxpar(head_oivis,'OI_REVN')
 date_obs = fxpar(head_oivis,'DATE-OBS')
 arrname  = fxpar(head_oivis,'ARRNAME')
 insname  = fxpar(head_oivis,'INSNAME')
 fxbread,unit,vis_target_id,'TARGET_ID'
 fxbread,unit,vis_time,'TIME'
 fxbread,unit,vis_mjd,'MJD'
 fxbread,unit,vis_int_time,'INT_TIME'
 fxbread,unit,visamp,'VISAMP'
 fxbread,unit,visamperr,'VISAMPERR'
 fxbread,unit,visphi,'VISPHI'
 fxbread,unit,visphierr,'VISPHIERR'
 fxbread,unit,vis_ucoord,'UCOORD'
 fxbread,unit,vis_vcoord,'VCOORD'
 fxbread,unit,vis_staindex,'STA_INDEX'
 fxbread,unit,vis_flag,'FLAG'

 nwave=n_elements(vis_flag)/nvis  ; WAVE implicitly defined 

 define_oivis,oivis,nwave=nwave ; structure depends on nwave (bad?)
 oivis=replicate(oivis,nvis)
 oivis(*).extver = extver_vis
 oivis(*).oi_revn = oi_revn
 oivis(*).date_obs = date_obs
 oivis(*).arrname = arrname
 oivis(*).insname = insname
 oivis.target_id=vis_target_id
 oivis.time=vis_time
 oivis.mjd =vis_mjd
 oivis.int_time=vis_int_time
 oivis.ucoord=vis_ucoord
 oivis.vcoord=vis_vcoord
 oivis.sta_index=vis_staindex
 for  i=0,nvis-1 do begin
  if ( nwave gt 1) then begin ; 2-d array (nwave >1)
   oivis(i).visamp=ptr_new(visamp(*,i))
   oivis(i).visamperr=ptr_new(visamperr(*,i))
   oivis(i).visphi=ptr_new(visphi(*,i))
   oivis(i).visphierr=ptr_new(visphierr(*,i))
   oivis(i).flag=ptr_new(vis_flag(*,i))
  endif else begin
   oivis(i).visamp=ptr_new(visamp(i))
   oivis(i).visamperr=ptr_new(visamperr(i))
   oivis(i).visphi=ptr_new(visphi(i))
   oivis(i).visphierr=ptr_new(visphierr(i))
   oivis(i).flag=ptr_new(vis_flag(i))
  endelse
 endfor ; make sure looping over right variable...

; If previous array exists, then append; otherwise create it
  if (n_elements(oivis_full) eq 0) then oivis_full=oivis $
   else oivis_full=concat_oitable(oivis_full,oivis)

endif ; OI_VIS

if (extname eq 'OI_VIS2') then begin   ; Read in OI_VIS2
 extver_vis2=extver_vis2+1
 head_oivis2=header

 nvis2=fxpar(head_oivis2,'NAXIS2')
 oi_revn =  fxpar(head_oivis2,'OI_REVN')
 date_obs = fxpar(head_oivis2,'DATE-OBS')
 arrname  = fxpar(head_oivis2,'ARRNAME')
 insname  = fxpar(head_oivis2,'INSNAME')
 fxbread,unit,vis2_target_id,'TARGET_ID'
 fxbread,unit,vis2_time,'TIME'
 fxbread,unit,vis2_mjd,'MJD'
 fxbread,unit,vis2_int_time,'INT_TIME'
 fxbread,unit,vis2data,'VIS2DATA'
 fxbread,unit,vis2err,'VIS2ERR'
 fxbread,unit,vis2_ucoord,'UCOORD'
 fxbread,unit,vis2_vcoord,'VCOORD'
 fxbread,unit,vis2_staindex,'STA_INDEX'
 fxbread,unit,vis2_flag,'FLAG'

 nwave=n_elements(vis2_flag)/nvis2  ; mplicitly defined (consistent?)

 define_oivis2,oivis2,nwave=nwave ; structure depends on nwave (bad?)
 oivis2=replicate(oivis2,nvis2)
 oivis2(*).extver = extver_vis
 oivis2(*).oi_revn = oi_revn
 oivis2(*).date_obs = date_obs
 oivis2(*).arrname = arrname
 oivis2(*).insname = insname
 oivis2.target_id=vis2_target_id
 oivis2.time=vis2_time
 oivis2.mjd =vis2_mjd
 oivis2.int_time=vis2_int_time
 oivis2.ucoord=vis2_ucoord
 oivis2.vcoord=vis2_vcoord
 oivis2.sta_index=vis2_staindex
for  i=0,nvis2-1 do begin
 if (nwave gt 1) then begin
  oivis2(i).vis2data=ptr_new(vis2data(*,i))
  oivis2(i).vis2err=ptr_new(vis2err(*,i))
  oivis2(i).flag=ptr_new(vis2_flag(*,i))
 endif else begin
  oivis2(i).vis2data=ptr_new(vis2data(i))
  oivis2(i).vis2err=ptr_new(vis2err(i))
  oivis2(i).flag=ptr_new(vis2_flag(i))
 endelse
endfor

; If previous array exists, then append; otherwise create it
  if (n_elements(oivis2_full) eq 0) then oivis2_full=oivis2 $
   else oivis2_full=concat_oitable(oivis2_full,oivis2)

 endif ; OI_VIS2

if (extname eq 'OI_T3') then begin   ; Read in OI_T3
 extver_t3=extver_T3+1
 head_oit3=header

 nt3=fxpar(head_oit3,'NAXIS2')
 oi_revn =  fxpar(head_oit3,'OI_REVN')
 date_obs = fxpar(head_oit3,'DATE-OBS')
 arrname  = fxpar(head_oit3,'ARRNAME')
 insname  = fxpar(head_oit3,'INSNAME')
 fxbread,unit,t3_target_id,'TARGET_ID'
 fxbread,unit,t3_time,'TIME'
 fxbread,unit,t3_mjd,'MJD'
 fxbread,unit,t3_int_time,'INT_TIME'
 fxbread,unit,t3amp,'T3AMP'
 fxbread,unit,t3amperr,'T3AMPERR'
 fxbread,unit,t3phi,'T3PHI'
 fxbread,unit,t3phierr,'T3PHIERR'
 fxbread,unit,t3_u1coord,'U1COORD'
 fxbread,unit,t3_v1coord,'V1COORD'
 fxbread,unit,t3_u2coord,'U2COORD'
 fxbread,unit,t3_v2coord,'V2COORD'
 fxbread,unit,t3_staindex,'STA_INDEX'
 fxbread,unit,t3_flag,'FLAG'
 nwave=n_elements(t3_flag)/nt3
 
 define_oit3,oit3,nwave=nwave ; structure depends on nwave (bad?)
 oit3=replicate(oit3,nt3)
 oit3(*).extver  = extver_t3
 oit3(*).oi_revn= oi_revn
 oit3(*).date_obs = date_obs
 oit3(*).arrname = arrname
 oit3(*).insname = insname 
 oit3.target_id=t3_target_id
 oit3.time=t3_time
 oit3.mjd = t3_mjd
 oit3.int_time=t3_int_time
 oit3.u1coord=t3_u1coord
 oit3.v1coord=t3_v1coord
 oit3.u2coord=t3_u2coord
 oit3.v2coord=t3_v2coord
 oit3.sta_index=t3_staindex
 for  i=0,nt3-1 do begin
  if ( nwave gt 1) then begin
   oit3(i).t3amp=ptr_new(t3amp(*,i))
   oit3(i).t3amperr=ptr_new(t3amperr(*,i))
   oit3(i).t3phi=ptr_new(t3phi(*,i))
   oit3(i).t3phierr=ptr_new(t3phierr(*,i))
   oit3(i).flag=ptr_new(t3_flag(*,i))
  endif else begin
   oit3(i).t3amp=ptr_new(t3amp(i))
   oit3(i).t3amperr=ptr_new(t3amperr(i))
   oit3(i).t3phi=ptr_new(t3phi(i))
   oit3(i).t3phierr=ptr_new(t3phierr(i))
   oit3(i).flag=ptr_new(t3_flag(i))
  endelse
 endfor ; make sure looping over right variable...

; If previous array exists, then append; otherwise create it
  if (n_elements(oit3_full) eq 0) then oit3_full=oit3 $
   else oit3_full=concat_oitable(oit3_full,oit3)
endif ; OI_T3

 fxbclose,unit,errmsg=errmsg

exten=exten+1
fxbopen, unit,file,exten,header,errmsg=errmsg

;print,'errmsg: ',errmsg,' exten: ',exten

endwhile

; Now do check for required tables.

if (keyword_set(inventory) eq 1) then begin
if (extver_target eq 1 and extver_wavelength ge 1 and $
    (extver_vis2 ge 1 or extver_vis ge 1 or extver_t3 ge 1)) then $
  print,'This file Satisfies the requirements of the OI_DATA format'

print,"Inventory:"
print,"  OI_ARRAY:      ",extver_array
print,"  OI_TARGET:     ",extver_target
print,"  OI_WAVELENGTH: ",extver_wavelength
print,"  OI_VIS:        ",extver_vis
print,"  OI_VIS2:       ",extver_vis2
print,"  OI_T3:         ",extver_t3
print," Unknown Tables: ",unknown_table
endif

if (extver_array ge 1) then oiarray=oiarray_full
if (extver_target ge 1) then oitarget=oitarget_full
if (extver_wavelength ge 1) then oiwavelength=oiwavelength_full
if (extver_vis ge 1) then oivis=oivis_full
if (extver_vis2 ge 1) then oivis2=oivis2_full
if (extver_t3 ge 1) then oit3=oit3_full

end
