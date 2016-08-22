pro extract_t3data, file=file, oit3=oit3_input, $
 oiwavelength=oiwavelength_input, oitarget=oitarget_input,t3data  ,$
 status=status

;+
; NAME:
;     extract_uvdata
; PURPOSE:
;     To extract uvdata into 'usuable' IDL variables for plotting,
;     inspecting data, etc.  Can either be used to read directly from
;     FITS file or by passing the required OITABLES. Will also allow (someday)
;     selection parameters
;
; CALLING SEQUENCE (Examples):
;     extract_t3data, t3data, file=fits_file, $
;             eff_wave =[1.8,2.3], target_id =[1,2]
;           or
;     extract_t3data, t3data, oit3= oit3table,$
;             oiwavelength=oiwavetable, oitarget=oitargettable,$
;             insname=['IOTA_IONIC3']
;
; INPUTS:
;
;   KEYWORDS:

;   One must pass EITHER 
;     FILE      - Name of fits file to read. Must be in
;
;            or ALL
;   OIWAVELENGTH Contains information on the observing bandpasses.
;      OIT3     Contains visibility-squared measurements.
;      OITARGET  Contains information on the target names

;
;    -
;
; OPTIONAL INPUTS:
;
;     Various keywords can be used to select subsets of the data that
;     pass certain selection criteria.  This section will be expanded
;     and new capabilities are added.
;
;     /STATUS    Print out status messages as data is being Extracted
;
; OUTPUTS:
;
;      T3DATA    A structure which contains in a simple format all the
;                requested data.   Does not use pointers and so all
;                data can be easily searched and plotted in IDL.
;                See define_oit3.pro for detailed description of
;                this structure
;
; RESTRICTIONS:
;      None.
;
; PROCEDURE:
;      Call from command line or in a program.
;
;
;      The Data Exchange Standard for Optical (Visible/IR) Interferometry
;      is being maintained by the IAU Working group on Optical Interferometry
;      and the current version of the standard can be found at :
;            http://www.mrao.cam.ac.uk/~jsy1001/exchange/
;
; EXAMPLE:


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
;     v0.0 2003Sep26    J.D. Monnier    based on extract_vis2data
;     v0.1 2005May30    J.D. Monnier    Added target names.
;
;     To do: A.  Add Revision Checking (once multiple revisions exist)  
;	     B.  Add /gui or /interactive options to do fancy
;                selections.
;            C.  Better options for SORTING by TIME, baseline, etc.
;            D.  Better error checking
;  *** Please write monnier@umich.edu with bug reports. ***
;-

if (keyword_set(file) eq 1) then begin
  read_oidata, file, oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
endif else begin
  oiwavelength = oiwavelength_input
  oit3 =oit3_input
   oitarget=oitarget_input
endelse

; Add error checking.

; First approach [maybe too slow and memory hog for big data sets!]
;    make big long array.
;    apply cuts at the end using a WHERE statements.
;

if n_elements(oit3) eq 0 then begin
  return
endif

lastarrname = ' '
;init
define_t3data,t3data_unit
insnames = strtrim( oiwavelength.insname,2)
t3_insnames =strtrim(oit3.insname,2)

for i=0,n_elements(oit3)-1 do begin
  oiwave_index =  (where(insnames eq t3_insnames(i)))(0)
   nwave = oiwavelength(oiwave_index).nwave
       ; Not supposed to ever crash in the above line since
       ; oidata standard requires wavelength tables for each insnames

if (keyword_set(status) eq 1) then begin
   if (oit3(i).arrname ne lastarrname) then begin
       pe= strtrim(string(fix(100.*i/n_elements(oit3))),2)
       print, pe+'% Complete:  Beginning extraction of ARRNAME : ',oit3(i).arrname 
       lastarrname = oit3(i).arrname
   endif
endif

   t3_new = replicate(t3data_unit,nwave)
   
   t3_new(*).oi_revn = oit3(i).oi_revn
   t3_new(*).date_obs= oit3(i).date_obs
   t3_new(*).arrname = strtrim(oit3(i).arrname,2)
   t3_new(*).insname = strtrim(oit3(i).insname,2)
   t3_new.eff_wave   = *oiwavelength(oiwave_index).eff_wave
   t3_new.eff_band   = *oiwavelength(oiwave_index).eff_band
   tin=where( oitarget.target_id eq oit3(i).target_id,ct ) ; NO ERRORS!!
    if ct eq 0 then begin
        print,' Target ID not in OITARGET table!! 
	stop
        return
    endif
   t3_new(*).target     = oitarget(tin).target

   t3_new(*).target_id  = oit3(i).target_id
   t3_new(*).time       = oit3(i).time
   t3_new(*).mjd        = oit3(i).mjd
   t3_new(*).int_time   = oit3(i).int_time
   t3_new.t3amp      =*oit3(i).t3amp
   t3_new.t3amperr       =*oit3(i).t3amperr
   t3_new.t3phi      =*oit3(i).t3phi
   t3_new.t3phierr       =*oit3(i).t3phierr
   t3_new(*).u1coord     = oit3(i).u1coord
   t3_new(*).v1coord     = oit3(i).v1coord
   t3_new(*).u2coord     = oit3(i).u2coord
   t3_new(*).v2coord     = oit3(i).v2coord
   t3_new(*).sta_index   = oit3(i).sta_index
   t3_new.flag          =*oit3(i).flag
  
; do baseline geometry calc at the end..

   if (i eq 0) then t3data=t3_new else t3data=concat_oitable(t3data,t3_new)
   
endfor

; Geometry
   t3data.u3coord= -1.0*t3data.u1coord - t3data.u2coord
   t3data.v3coord= -1.0*t3data.v1coord - t3data.v2coord

   t3data.u1 = t3data.u1coord/t3data.eff_wave
   t3data.v1 = t3data.v1coord/t3data.eff_wave
   t3data.u2 = t3data.u2coord/t3data.eff_wave
   t3data.v2 = t3data.v2coord/t3data.eff_wave
   t3data.u3 = t3data.u3coord/t3data.eff_wave
   t3data.v3 = t3data.v3coord/t3data.eff_wave
   ri2at, t3data.v1coord, t3data.u1coord, amp1, theta1 
            ;not a bug (angle E of North)
 ri2at, t3data.v2coord, t3data.u2coord, amp2, theta2 ; 
 ri2at, t3data.v3coord, t3data.u3coord, amp3, theta3 ; 

   t3data.baseline1 = amp1
   t3data.pa1       = theta1
   t3data.baseline2 = amp2
   t3data.pa2       = theta2
   t3data.baseline3 = amp3
   t3data.pa3       = theta3


end


