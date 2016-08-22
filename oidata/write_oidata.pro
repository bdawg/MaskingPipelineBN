pro write_oidata, outfile,oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
;+
; NAME:
;     write_oidata
; PURPOSE:
;     To save FITS files in the OI-DATA format based on information in 
;      specific IDL structures.
;
; CALLING SEQUENCE:
;     write_oidata, fits_file, oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
;                       or
;     write_oidata, fits_file, oiarray,oitarget,oiwavelength,0,oivis2,0
;	       In the latter call, the oivis and oit3 tables are not written.
; INPUTS:
;     fits_file - Name of FITS file to saved  in OI-DATA format. 
;		  See Modification History to see what Revision #s are 
;		  supported and for URL of the OI-DATA format definitions.
;
;  Note: since multiple versions of most tables are allowed in this OI_DATA
;        format, this routine will group tables together based on the
;        ARRNAME, INSNAME, and ARRNAME/INSNAME for OI_ARRAY, OI_WAVELENGTH, 
;        and OI_VIS/OI_VIS2/OI_T3 respectively.  This also means that
;        data arrays which can have variable length (e.g., NWAVE) must be
;	 referred to by their pointers.  Hence, these structures contain
;	 *pointer* to arrays generally, not the arrays themselves.
;
;      OIARRAY   Contains information on the interferometer array geometry.
; 		 IDL array of OIARRAY structures containing all the information
;		 necessary for the OI_ARRAY Binary FITS Table.
;		 There is one entry for each row in the binary 
;		 table (one row for each telescope) and the common header 
;		 keywords are repeated in each row.  Different OIARRAY tables
;		 are distinguished using the ARRNAME (and EXTVER) keyword.
;      OITARGET	 Contains information on the observed targets. 
;                IDL array of OITARGET structures containing all the information
;                necessary for the OI_TARGET Binary FITS Table.
;                There is one entry for each row in the binary 
;                table (one row for each target) and the common header 
;                keywords are repeated in each row. Only one OITARGET table
;		 allowed.
;   OIWAVELENGTH Contains information on the observing bandpasses.
;                IDL array of OIWAVELENGTH structures containing all the 
;                information necessary for the OI_WAVELENGTH Binary FITS Table.
;                information in the OI_WAVELENGTH Binary FITS Table if present
;                in the OIDATA FITS file.  There is one entry for each instance
;                of the OIWAVELENGTH table and we note that EFF_WAVE/EFF_BAND
;                are pointers, not arrays themselves. Different OIWAVELENGTH
;                tables are differentiated also by different INSNAME keywords.
;      OIVIS	 Contains complex visibility measurements.
;		 IDL array of OIVIS structures containing all the 
;                information necessary for the OI_VIS Binary FITS Table.
;                There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row.  
;		 Note that the structures contain pointers to the data, not
;		 the data themselves. Also, different instances of each OIVIS
;		 table are differentiated using combination of ARRNAME/INSNAME
;		 and DATE-OBS.
;		 
;      OIVIS2    Contains visibility-squared measurements.
;                IDL array of OIVIS2 structures containing all the 
;                information necessary for the OI_VIS2 Binary FITS Table.
;                There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row.
;		 Note that the structures contain pointers to the data, not
;                the data themselves. Also, different instances of each OIVIS2
;                table are differentiated using combination of ARRNAME/INSNAME
;		 and DATE-OBS.

;      OIT3      Contains measurements of the bispectrum (triple amp and phase).
;                IDL array of OIT3 structures containing all the 
;                information necessary for the OI_T3 Binary FITS Table.
;                There is one entry for each row in 
;                the binary table (one row for each measurement) and the common 
;                header keywords are repeated in each row.
;                Note that the structures contain pointers to the data, not
;                the data themselves. Also, different instances of each OIT3
;                table are differentiated using combination of ARRNAME/INSNAME
;		 and DATE-OBS.
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;     None
; 
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      There are currently 6 different binary tables defined in the OIDATA
;      format.  All tables are not required to be present to be a valid OIDATA 
;      FITS file.  WRITE_OIDATA will write a FITS file containing the
;      any or all of 6 allowed tables.  The data for each table is read from the
;      IDL structure arrays.  The structures are defined in procedures
;      DEFINE_OI*.PRO elsewhere in the IOI Library. 
;
;      One potentially confusing aspect of the OI-DATA standard is that
;      multiple instances of OIARRAY, OI-WAVELENGTH, OI-VIS2, OI-VIS, and OI-T3
;      tables (all but OI-TARGET) are allowed, to encourage merging of data 
;      from different times and different interferometers. In order to
;      use structures, one must use Pointers to variable-length arrays, rather
;      than the arrays themselves.  Hence, all the arrays of NWAVE length
;      (e.g., the data and detector characteristics) are referenced by pointers.
;
;      The Data Exchange Standard for Optical (Visible/IR) Interferometry
;      is being maintained by the IAU Working group on Optical Interferometry
;      and the current version of the standard can be found at :
;            http://www.mrao.cam.ac.uk/~jsy1001/exchange/

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
;
;	IDL> WRITE_OIDATA, 'testdata1.fits', oiarray,oitarget,oiwavelength,$
;		oivis, oivis2,oit3
;	IDL> print,oivis2(0).time,*oivis2(0).vis2data,*oivis2(0).vis2err
; 		==>        82810.000      0.67700000     0.064000000
;       IDL> print,string(*oivis(0).flag)
;		==> 	   F
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
;	   included support of multiple instances of OI-Tables.
;     v0.2 2003Feb13    J.D. monnier    Updated to reflect oi-TARGET changes
;     v0.3 2003Sep26    J.D. Monnier    Small bug fix in datatype
;
;     To do:   Error checking if all columns headings and keywords are not used
;              (should fill in default dummy values)
;
;  *** Please write monnier@umich.edu with bug reports. ***
;-

; Modified by BRN Sept 2013- change loop type to LONG

;---------------------------------------------------------------------------
; Check to see if/which IDL structures are being passed. If not structures,
; then we assume this is a dummy variable
oi_types=[size(oiarray,/type), size(oitarget,/type), size(oiwavelength,/type),$
          size(oivis,/type)  , size(oivis2,/type)  , size(oit3,/type)]
; Type 8 indicates a structure

if (oi_types(1) ne 8 or oi_types(2) ne 8 $
    or (oi_types(3) ne 8 and oi_types(4) ne 8 and oi_types(5) ne 8)) then begin
 print,' All required binary tables for OI_DATA format are not present 
 print,' We will write the file anyway for now.  '
endif else begin
 
;print,' The essential binary tables for OI_DATA format are present  '
endelse

; NOW WRITE IT!
fxhmake, newhead_hdu, /init,/extend,fix(0)
fxaddpar, newhead_hdu, "COMMENT", $
 "  FITS (Flexible Image Transport System) format defined in Astronomy and"
fxaddpar, newhead_hdu, "COMMENT", $
 "  Astrophysics Supplement Series v44/p363, v44/p371, v73/p359, v73/p365."
fxaddpar, newhead_hdu, "COMMENT", $
 "  Contact the NASA Science Office of Standards and Technology for the   "
fxaddpar, newhead_hdu, "COMMENT", $
 "  FITS Definition document #100 and other FITS information.             "
fxaddpar, newhead_hdu, "COMMENT", $
 "  Binary Extensions conform to OI-DATA standard for exchange of "
fxaddpar, newhead_hdu, "COMMENT", $
 "  optical interferometry data currently described at"
fxaddpar, newhead_hdu, "COMMENT", $
 "  http://www.mrao.cam.ac.uk/~jsy1001/exchange/."
fxaddpar, newhead_hdu, "COMMENT", $
 "  This file was created by WRITE_OIDATA.PRO by IDL (library can be"
fxaddpar, newhead_hdu, "COMMENT", $
 "  obtained from www.astro.lsa.umich.edu/~monnier). "

; Create fits file
fxwrite,outfile,newhead_hdu

; OIARRAY
if (oi_types(0) eq 8) then begin

; We must allow for mutiple ARRAY tables..
 uniq_arrname = uniq(oiarray.arrname,sort(oiarray.arrname) )
 

for i_a=0,n_elements(uniq_arrname)-1 do begin
  uniq_index = where(oiarray.arrname eq oiarray(uniq_arrname(i_a)).arrname)
  oiarray1=oiarray(uniq_index)

  ntel=n_elements(oiarray1)
  fxbhmake, newhead_oiarray, ntel, 'OI_ARRAY','name of this binary table extension',EXTVER=i_a+1 ;extver >= 1
  fxbaddcol,acol,newhead_oiarray, "                ", 'TEL_NAME'
  fxbaddcol,bcol,newhead_oiarray, "                ", 'STA_NAME'
  fxbaddcol,ccol,newhead_oiarray, 1, 'STA_INDEX'
  fxbaddcol,dcol,newhead_oiarray, 1.0, 'DIAMETER',TUNIT='m'
  fxbaddcol,ecol,newhead_oiarray, dblarr(3) , 'STAXYZ',tunit='m'

  fxaddpar, newhead_oiarray, "OI_REVN", oiarray1(0).oi_revn,  $
     "Revision number of the table definition"
  fxaddpar, newhead_oiarray, "ARRNAME", oiarray1(0).arrname, $
     "Array name"
  fxaddpar, newhead_oiarray, "FRAME", oiarray1(0).frame, $
     "Coordinate frame"
  fxaddpar, newhead_oiarray, "ARRAYX", oiarray1(0).arrayx, $
     "[m] Array center x coordinate"
  fxaddpar, newhead_oiarray, "ARRAYY", oiarray1(0).arrayy, $
     "[m] Array center y coordinate"
  fxaddpar, newhead_oiarray, "ARRAYZ", oiarray1(0).arrayz, $
     "[m] Array center z coordinate"
  fxbcreate,unit,outfile,newhead_oiarray
  fxbwritm, unit, ["TEL_NAME","STA_NAME","STA_INDEX","DIAMETER","STAXYZ"], oiarray1.tel_name, $
     oiarray1.sta_name, oiarray1.sta_index,oiarray1.diameter,oiarray1.staxyz
  fxbfinish,unit
endfor ; loop over all OI_ARRAY tables.
endif 

; OITARGET
if (oi_types(1) eq 8) then begin
; There can only be one OI_TARGET table.
  ntargets=n_elements(oitarget)
  fxbhmake, newhead_oitarget, ntargets,$
    'OI_TARGET','name of this binary table extension'
  fxbaddcol,acol,newhead_oitarget, 1, 'TARGET_ID'
  fxbaddcol,bcol,newhead_oitarget, "                ", 'TARGET'
  fxbaddcol,ccol,newhead_oitarget, 1d0, 'RAEP0',tunit="deg"
  fxbaddcol,dcol,newhead_oitarget, 1d0, 'DECEP0',TUNIT='deg'
  fxbaddcol,ecol,newhead_oitarget, 1.0 , 'EQUINOX',tunit='yr'
  fxbaddcol,hcol,newhead_oitarget, 1d0, 'RA_ERR',tunit='deg'
  fxbaddcol,icol,newhead_oitarget, 1d0, 'DEC_ERR',TUNIT='deg'
  fxbaddcol,jcol,newhead_oitarget, 1d0 , 'SYSVEL',tunit='m/s'
  fxbaddcol,kcol,newhead_oitarget, "        ", 'VELTYP'
  fxbaddcol,lcol,newhead_oitarget, "        ", 'VELDEF'
  fxbaddcol,mcol,newhead_oitarget, 1d0, 'PMRA',tunit='deg/yr'
  fxbaddcol,ocol,newhead_oitarget, 1d0 , 'PMDEC',tunit='deg/yr'
  fxbaddcol,ncol,newhead_oitarget, 1d0, 'PMRA_ERR',TUNIT='deg/yr'
  fxbaddcol,pcol,newhead_oitarget, 1d0 , 'PMDEC_ERR',tunit='deg/yr'
  fxbaddcol,qcol,newhead_oitarget, 1.0, "PARALLAX", tunit='deg'
  fxbaddcol,rcol,newhead_oitarget, 1.0, 'PARA_ERR',tunit='deg'
  fxbaddcol,scol,newhead_oitarget, "                ", "SPECTYP"
 
  fxaddpar, newhead_oitarget, "OI_REVN", oitarget(0).oi_revn,  $
     "Revision number of the table definition"
  fxbcreate,unit,outfile,newhead_oitarget
  fxbwritm, unit, ["TARGET_ID","TARGET","RAEP0","DECEP0","EQUINOX",$
      "RA_ERR","DEC_ERR","SYSVEL","VELTYP","VELDEF","PMRA","PMDEC",$
      "PMRA_ERR","PMDEC_ERR","PARALLAX","PARA_ERR","SPECTYP"],$
    oitarget.target_id, oitarget.target,oitarget.raep0,oitarget.decep0,$
    oitarget.equinox, oitarget.ra_err,$
    oitarget.dec_err, oitarget.sysvel, oitarget.veltyp,oitarget.veldef,$
    oitarget.pmra,oitarget.pmdec,oitarget.pmra_err,oitarget.pmdec_err,$
    oitarget.parallax,oitarget.para_err,oitarget.spectyp

  fxbfinish,unit
endif

; OI_WAVELENGTH
if (oi_types(2) eq 8) then begin

; We must allow for mutiple OI_WAVELENGTH tables..
 uniq_insname = uniq(oiwavelength.insname,sort(oiwavelength.insname) )

for i_a=0,n_elements(uniq_insname)-1 do begin
  uniq_index = where(oiwavelength.insname eq oiwavelength(uniq_insname(i_a)).insname)
  oiwavelength1=oiwavelength(uniq_index)

  nwave=oiwavelength1(0).nwave
  fxbhmake, newhead_oiwavelength, nwave,$
    'OI_WAVELENGTH','name of this binary table extension',EXTVER=i_a+1
  fxbaddcol,acol,newhead_oiwavelength, 1.0, 'EFF_WAVE', tunit='m'
  fxbaddcol,bcol,newhead_oiwavelength, 1.0, 'EFF_BAND', tunit='m'
  
  fxaddpar, newhead_oiwavelength, "OI_REVN", oiwavelength1(0).oi_revn,  $
     "Revision number of the table definition"
  fxaddpar, newhead_oiwavelength, "INSNAME", oiwavelength1(0).insname,  $
     "Name of detector, for cross-referencing"
  fxbcreate,unit,outfile,newhead_oiwavelength
  fxbwritm, unit, ["EFF_WAVE","EFF_BAND"],$
    float(*oiwavelength1.eff_wave),float(*oiwavelength1.eff_band)

  fxbfinish,unit
 endfor ; Loop through all instances of the OI_WAVELENGTH array.
endif ;OI_WAVELENGTH

; OIVIS
if (oi_types(3) eq 8) then begin

; We must allow for mutiple OIVIS tables..
; Find unique arrname, visname, date-obs
 uniq_arrname = uniq(oivis.arrname,sort(oivis.arrname) )
 uniq_insname = uniq(oivis.insname,sort(oivis.insname) )
 uniq_dateobs = uniq(oivis.date_obs,sort(oivis.date_obs) )

extver=1
for i_a=0,n_elements(uniq_arrname)-1 do begin
 for j_a=0,n_elements(uniq_insname)-1 do begin
   for k_a=0,n_elements(uniq_dateobs)-1 do begin

  uniq_index = where(oivis.arrname eq oivis(uniq_arrname(i_a)).arrname  and $
                     oivis.insname eq oivis(uniq_insname(j_a)).insname  and $
                     oivis.date_obs eq oivis(uniq_dateobs(k_a)).date_obs,ct_left)
  if (ct_left ne 0) then begin

  oivis1=oivis(uniq_index)

  nvis=n_elements(oivis1)
  nwave = oivis1(0).nwave
  fxbhmake, newhead_oivis, nvis, $
    'OI_VIS','name of this binary table extension',extver=extver
  fxbaddcol,acol,newhead_oivis, 1, 'TARGET_ID'
  fxbaddcol,bcol,newhead_oivis, 1d0, 'TIME',tunit='s'
  fxbaddcol,ccol,newhead_oivis, 1d0, 'MJD',tunit='day'
  fxbaddcol,dcol,newhead_oivis, 1d0, 'INT_TIME',tunit='s'
  fxbaddcol,ecol,newhead_oivis, dindgen(nwave), 'VISAMP',"Visibility Amplitude"
  fxbaddcol,fcol,newhead_oivis, dindgen(nwave), 'VISAMPERR'
  fxbaddcol,gcol,newhead_oivis, dindgen(nwave), 'VISPHI'
  fxbaddcol,hcol,newhead_oivis, dindgen(nwave), 'VISPHIERR'
  fxbaddcol,icol,newhead_oivis, 1d0, 'UCOORD', tunit='m'
  fxbaddcol,jcol,newhead_oivis, 1d0, 'VCOORD',tunit='m'
  fxbaddcol,jcol,newhead_oivis, [1,2], 'STA_INDEX'
  fxbaddcol,lcol,newhead_oivis, bytarr(nwave), 'FLAG', /logical
 
  fxaddpar, newhead_oivis, "OI_REVN", oivis1(0).oi_revn,  $
     "Revision number of the table definition"
  fxaddpar, newhead_oivis, "DATE-OBS", oivis1(0).date_obs,  $
     "UTC start date of observations"
  fxaddpar, newhead_oivis, "ARRNAME", oivis1(0).arrname, $
     "(optional) Identifies corresponding OI_ARRAY"
  fxaddpar, newhead_oivis, "INSNAME", oivis1(0).insname, $
     "Identifies corresponding OI_WAVELENGTH table"

  fxbcreate,unit,outfile,newhead_oivis
 ;Setup DATA ARRAYS in order to pass data properly/efficiently
  visamp = dblarr(nwave,nvis)
  visamperr= dblarr(nwave,nvis)
  visphi   = dblarr(nwave,nvis)
  visphierr= dblarr(nwave,nvis)
  visflag  = bytarr(nwave,nvis)  
  for i=0,nvis-1 do begin
   if (nwave gt 1) then begin 
    visamp(*,i) = *oivis1(i).visamp
    visamperr(*,i) = *oivis1(i).visamperr
    visphi(*,i) = *oivis1(i).visphi
    visphierr(*,i) = *oivis1(i).visphierr
    visflag(*,i)   = *oivis1(i).flag
   endif else begin
    visamp(i) = *oivis1(i).visamp
    visamperr(i) = *oivis1(i).visamperr
    visphi(i) = *oivis1(i).visphi
    visphierr(i) = *oivis1(i).visphierr
    visflag(i)   = *oivis1(i).flag
   endelse
  endfor
        
  fxbwritm, unit, ["TARGET_ID","TIME","MJD","INT_TIME","VISAMP","VISAMPERR",$
      "VISPHI","VISPHIERR","UCOORD",$
      "VCOORD","STA_INDEX","FLAG"],$
      oivis1.target_id,oivis1.time,oivis1.mjd,oivis1.int_time,visamp,visamperr,$
      visphi,visphierr,oivis1.ucoord,oivis1.vcoord,oivis1.sta_index,visflag
  fxbfinish,unit
  extver=extver+1
endif
endfor
endfor
endfor
endif

; OIVIS2
if (oi_types(4) eq 8) then begin

; We must allow for mutiple OIVIS2 tables..
; Find unique arrname, visname, date-obs
 uniq_arrname = uniq(oivis2.arrname,sort(oivis2.arrname) )
 uniq_insname = uniq(oivis2.insname,sort(oivis2.insname) )
 uniq_dateobs = uniq(oivis2.date_obs,sort(oivis2.date_obs) )

extver=1
for i_a=0,n_elements(uniq_arrname)-1 do begin
 for j_a=0,n_elements(uniq_insname)-1 do begin
   for k_a=0,n_elements(uniq_dateobs)-1 do begin

  uniq_index = where(oivis2.arrname eq oivis2(uniq_arrname(i_a)).arrname  and $
                     oivis2.insname eq oivis2(uniq_insname(j_a)).insname  and $
                     oivis2.date_obs eq oivis2(uniq_dateobs(k_a)).date_obs,ct_left)
  if (ct_left ne 0) then begin

  oivis1=oivis2(uniq_index)


  nvis2=n_elements(oivis1)
  nwave = oivis1(0).nwave
  fxbhmake, newhead_oivis2, nvis2, $
    'OI_VIS2','name of this binary table extension',extver=extver
  fxbaddcol,acol,newhead_oivis2, 1, 'TARGET_ID'
  fxbaddcol,bcol,newhead_oivis2, 1d0, 'TIME',tunit='s'
  fxbaddcol,ccol,newhead_oivis2, 1d0, 'MJD',tunit='day'
  fxbaddcol,dcol,newhead_oivis2, 1d0, 'INT_TIME',tunit='s'
  fxbaddcol,ecol,newhead_oivis2, dindgen(nwave), 'VIS2DATA'
  fxbaddcol,fcol,newhead_oivis2, dindgen(nwave), 'VIS2ERR'
  fxbaddcol,gcol,newhead_oivis2, 1d0, 'UCOORD', tunit='m'
  fxbaddcol,hcol,newhead_oivis2, 1d0, 'VCOORD',tunit='m'
  fxbaddcol,icol,newhead_oivis2, [1,2], 'STA_INDEX'
  fxbaddcol,jcol,newhead_oivis2, bytarr(nwave), 'FLAG', /logical
 
 fxaddpar, newhead_oivis2, "OI_REVN", oivis1(0).oi_revn,  $
     "Revision number of the table definition"
  fxaddpar, newhead_oivis2, "DATE-OBS", oivis1(0).date_obs,  $
     "UTC start date of observations"
  fxaddpar, newhead_oivis2, "ARRNAME", oivis1(0).arrname, $
     "(optional) Identifies corresponding OI_ARRAY"
  fxaddpar, newhead_oivis2, "INSNAME", oivis1(0).insname, $
     "Identifies corresponding OI_WAVELENGTH table"

 ;Setup DATA ARRAYS in order to pass data properly/efficiently
  vis2data = dblarr(nwave,nvis2)
  vis2err  = dblarr(nwave,nvis2)
  vis2flag = bytarr(nwave,nvis2)
  for i=0,nvis2-1 do begin
   if (nwave gt 1) then begin
    vis2data(*,i) = *oivis1(i).vis2data
    vis2err(*,i)  = *oivis1(i).vis2err
    vis2flag(*,i) = *oivis1(i).flag
   endif else begin
    vis2data(i) = *oivis1(i).vis2data
    vis2err(i)  = *oivis1(i).vis2err
    vis2flag(i) = *oivis1(i).flag
   endelse
  endfor
  fxbcreate,unit,outfile,newhead_oivis2
  fxbwritm, unit, ["TARGET_ID","TIME","MJD","INT_TIME","VIS2DATA","VIS2ERR","UCOORD",$
      "VCOORD","STA_INDEX","FLAG"],$
      oivis1.target_id,oivis1.time,oivis1.mjd,oivis1.int_time,vis2data,vis2err,$
      oivis1.ucoord,oivis1.vcoord,oivis1.sta_index,vis2flag
  fxbfinish,unit
extver=extver+1
endif
endfor
endfor
endfor
endif ; DONE WRITING OI_VIS2


; OIT3
if (oi_types(5) eq 8) then begin

; We must allow for mutiple OIVIS2 tables..
; Find unique arrname, visname, date-obs
 uniq_arrname = uniq(oit3.arrname,sort(oit3.arrname) )
 uniq_insname = uniq(oit3.insname,sort(oit3.insname) )
 uniq_dateobs = uniq(oit3.date_obs,sort(oit3.date_obs) )

extver=1
for i_a=0,n_elements(uniq_arrname)-1 do begin
 for j_a=0,n_elements(uniq_insname)-1 do begin
   for k_a=0,n_elements(uniq_dateobs)-1 do begin

  uniq_index = where(oit3.arrname eq oit3(uniq_arrname(i_a)).arrname  and $
                     oit3.insname eq oit3(uniq_insname(j_a)).insname  and $
                     oit3.date_obs eq oit3(uniq_dateobs(k_a)).date_obs,ct_left)
  if (ct_left ne 0) then begin

  oita=oit3(uniq_index)


  nt3=n_elements(oita)
  nwave = oita(0).nwave
  fxbhmake, newhead_oit3, nt3, $
    'OI_T3','name of this binary table extension',extver=extver
  fxbaddcol,acol,newhead_oit3, 1, 'TARGET_ID'
  fxbaddcol,bcol,newhead_oit3, 1d0, 'TIME',tunit='s'
  fxbaddcol,ccol,newhead_oit3, 1d0, 'MJD',tunit='day'
  fxbaddcol,dcol,newhead_oit3, 1d0, 'INT_TIME',tunit='s'
  fxbaddcol,ecol,newhead_oit3, dindgen(nwave), 'T3AMP'
  fxbaddcol,fcol,newhead_oit3, dindgen(nwave), 'T3AMPERR'
  fxbaddcol,gcol,newhead_oit3, dindgen(nwave), 'T3PHI'
  fxbaddcol,hcol,newhead_oit3, dindgen(nwave), 'T3PHIERR'
  fxbaddcol,icol,newhead_oit3, 1d0, 'U1COORD', tunit='m'
  fxbaddcol,jcol,newhead_oit3, 1d0, 'V1COORD',tunit='m'
  fxbaddcol,kcol,newhead_oit3, 1d0, 'U2COORD', tunit='m'
  fxbaddcol,lcol,newhead_oit3, 1d0, 'V2COORD',tunit='m'
  fxbaddcol,mcol,newhead_oit3, [1,2,3], 'STA_INDEX'
  fxbaddcol,ncol,newhead_oit3, bytarr(nwave), 'FLAG', /logical
  fxaddpar, newhead_oit3, "OI_REVN", oita(0).oi_revn,  $
     "Revision number of the table definition"
  fxaddpar, newhead_oit3, "DATE-OBS", oita(0).date_obs,  $
     "UTC start date of observations"
  fxaddpar, newhead_oit3, "ARRNAME", oita(0).arrname, $
     "(optional) Identifies corresponding OI_ARRAY"
  fxaddpar, newhead_oit3, "INSNAME", oita(0).insname, $
     "Identifies corresponding OI_WAVELENGTH table"
 ;Setup DATA ARRAYS in order to pass data properly/efficiently
  t3amp = dblarr(nwave,nt3)
  t3amperr= dblarr(nwave,nt3)
  t3phi   = dblarr(nwave,nt3)
  t3phierr= dblarr(nwave,nt3)
  t3flag  = bytarr(nwave,nt3)  
  for i=0L,nt3-1 do begin
   if (nwave gt 1) then begin 
    t3amp(*,i) = *oita(i).t3amp
    t3amperr(*,i) = *oita(i).t3amperr
    t3phi(*,i) = *oita(i).t3phi
    t3phierr(*,i) = *oita(i).t3phierr
    t3flag(*,i)   = *oita(i).flag
   endif else begin
    t3amp(i) = *oita(i).t3amp
    t3amperr(i) = *oita(i).t3amperr
    t3phi(i) = *oita(i).t3phi
    t3phierr(i) = *oita(i).t3phierr
    t3flag(i)   = *oita(i).flag
   endelse
  endfor
        
  fxbcreate,unit,outfile,newhead_oit3
  fxbwritm, unit, ["TARGET_ID","TIME","MJD","INT_TIME","T3AMP","T3AMPERR",$
     "T3PHI","T3PHIERR", "U1COORD","V1COORD", "U2COORD","V2COORD",$
      "STA_INDEX","FLAG"],$
      oita.target_id,oita.time,oita.mjd,oita.int_time,t3amp,t3amperr,$
      t3phi,t3phierr,oita.u1coord,oita.v1coord,$
      oita.u2coord,oita.v2coord,oita.sta_index,t3flag
  fxbfinish,unit
 extver=extver+1
endif
endfor
endfor
endfor

endif ; DONE WRITING OI_T3


end
