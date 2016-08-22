pro extract_vis2data, file=file, oivis2=oivis2_input, $
 oiwavelength=oiwavelength_input, oitarget=oitarget_input,vis2data  ,$
 status=status

;+
; NAME:
;     extract_uvdata
; PURPOSE:
;     To extract uvdata into 'usuable' IDL variables for plotting,
;     inspecting data, etc.  Can either be used to read directly from
;     FITS file or by passing the required OITABLES. Will also allow 
;     selection parameters
;
; CALLING SEQUENCE (Examples):
;     extract_vis2data, vis2data, file=fits_file, $
;             eff_wave =[1.8,2.3], target_id =[1,2]
;           or
;     extract_vis2data, vis2data, oivis2= oivis2table,$
;             oiwavelength=oiwavetable, oitarget=oitargettable,$
;              insname=['IOTA_IONIC3']
;
; INPUTS:
;
;   KEYWORDS:

;   One must pass EITHER 
;     FILE      - Name of fits file to read. Must be in
;
;            or All
;   OIWAVELENGTH Contains information on the observing bandpasses.
;      OIVIS2    Contains visibility-squared measurements.
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
;      VIS2DATA  A structure which contains in a simple format all the
;                requested data.   Does not use pointers and so all
;                data can be easily searched and plotted in IDL.
;                See define_vis2data.pro for detailed description of
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
;     v0.0 2003Jul04    J.D. Monnier    Initiated
;     v0.1 2005May30    J.D. Monnier	Added target names.
;     v0.2 2006Oct22    J.D. Monnier	Added sfu
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
  oivis2 =oivis2_input
  oitarget=oitarget_input
endelse

; Add error checking.

; First approach [maybe too slow and memory hog for big data sets!]
;    make big long array.
;    apply cuts at the end using a WHERE statements.
;

lastarrname = ' '
;init
define_vis2data,vis2data_unit
insnames = strtrim( oiwavelength.insname,2)
vis2_insnames =strtrim(oivis2.insname,2)

for i=0,n_elements(oivis2)-1 do begin
  oiwave_index =  (where(insnames eq vis2_insnames(i)))(0)
   nwave = oiwavelength(oiwave_index).nwave
       ; Not supposed to ever crash in the above line since
       ; oidata standard requires wavelength tables for each insnames

if (keyword_set(status) eq 1) then begin
   if (oivis2(i).arrname ne lastarrname) then begin
       pe= strtrim(string(fix(100.*i/n_elements(oivis2))),2)
       print, pe+'% Complete:  Beginning extraction of ARRNAME : ',oivis2(i).arrname 
       lastarrname = oivis2(i).arrname
   endif
endif

   vis2_new = replicate(vis2data_unit,nwave)
   
   vis2_new(*).oi_revn = oivis2(i).oi_revn
   vis2_new(*).date_obs= oivis2(i).date_obs
   vis2_new(*).arrname = strtrim(oivis2(i).arrname,2)
   vis2_new(*).insname = strtrim(oivis2(i).insname,2)
   vis2_new.eff_wave   = *oiwavelength(oiwave_index).eff_wave
   vis2_new.eff_band   = *oiwavelength(oiwave_index).eff_band
   vis2_new(*).target_id  = oivis2(i).target_id
   tin=where( oitarget.target_id eq oivis2(i).target_id,ct ) ; NO ERRORS!!
    if ct eq 0 then begin
        print,' Target ID not in OITARGET table!! 
	stop
        return
    endif
   vis2_new(*).target     = oitarget(tin).target
   vis2_new(*).time       = oivis2(i).time
   vis2_new(*).mjd        = oivis2(i).mjd
   vis2_new(*).int_time   = oivis2(i).int_time
   vis2_new.vis2data      =*oivis2(i).vis2data
   vis2_new.vis2err       =*oivis2(i).vis2err
   vis2_new(*).ucoord     = oivis2(i).ucoord
   vis2_new(*).vcoord     = oivis2(i).vcoord
   vis2_new(*).sta_index  = oivis2(i).sta_index
   vis2_new.flag          =*oivis2(i).flag
  
; do baseline geometry calc at the end..

   if (i eq 0) then vis2data=vis2_new else vis2data=concat_oitable(vis2data,vis2_new)
   
endfor

; Geometry
   vis2data.u = vis2data.ucoord/vis2data.eff_wave
   vis2data.v = vis2data.vcoord/vis2data.eff_wave
   ri2at, vis2data.vcoord, vis2data.ucoord, amp, theta ; not a bug (angle E of No)
   vis2data.baseline = amp
   vis2data.pa       = theta
   vis2data.sfu = vis2data.baseline/vis2data.eff_wave

end


