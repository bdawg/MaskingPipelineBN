pro define_t3data,t3data_unit
;+
; NAME:
;     define_t3data
; PURPOSE:
;     Define the T3DATA structure which will contain the data in 
;     the OIDATA format in a easily accessible format for by IDL programs.
;
; CALLING SEQUENCE:
;     define_t3data, t3data_unit
;
; INPUTS:
;     None
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;      T3DATA_UNIT  
;		 Will contain measurements on the interferometer complex
;		 visibility. This structure can be read by WRITE_OIDATA
;		 to write standard OI-DATA FITS files.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      Used by extract_t3data.
;
;      The Data Exchange Standard for Optical (Visible/IR) Interferometry
;      is being maintained by the IAU Working group on Optical Interferometry
;      and the current version of the standard can be found at :
;            http://www.mrao.cam.ac.uk/~jsy1001/exchange/

; EXAMPLE:
;       This procedure is not meant to be called from the commandline.

; PROCEDURES USED:
;	This routine is part of the IDL Optical Interferometry (IOI) Library.
;        (more information at www.astro.lsa.umich.edu/~monnier)
;       The routines of the IOI Library generically 
;       require an up-to-date Astrolib library in order to read and write binary
;       FITS tables, among other things. The IDL Astronomy User's Library can
;       be found at http://idlastro.gsfc.nasa.gov/homepage.html.
;	
;
; MODIFICATION HISTORY:
;     v0.0 2003Sep26    J.D. Monnier    based on define vis2data
;     v0.1 2005May30	JDM		added target
;
;     To do: A.  Add Revision Keyword (once multiple revisions exist)  
;-

; The variables below are standard names are defined in the OIDATA standard,
; unless accompanied by IDL comment explaining defition.
;
; This data structure contains extra information that may not be necessary
; to replicate for EVERY datum, and thus may be cumbersome to use on 
; large data sets.  See extract_vis2data for more information.
;
t3data_unit = { $
  oi_revn: 1	,$ 
  date_obs: " "	,$ 
  arrname:  " "	,$
  insname:  " " ,$
  eff_wave: 1d0, $
  eff_band: 1d0, $
  target_id: 0	,$
  target: ' '   ,$
  time: 0d0	,$
  mjd:  0d0	,$
  int_time: 0d0	,$
  t3amp: 0d0 	,$ 
  t3amperr: 0d0, $
  t3phi:   0d0 ,$
  t3phierr: 0d0, $
  u1coord : 0d0	,$  
  v1coord : 0d0	,$
  u1      : 0d0  ,$   u in inverse radians !!
  v1      : 0d0  ,$   v in inverse radians !!
  baseline1: 0d0 ,$   Projected Baseline Length in Meters
  pa1      : 0d0 ,$   Position Angle of Baseline (Degrees E of N)
  u2coord : 0d0 ,$
  v2coord : 0d0 ,$
  u2      : 0d0  ,$   u in inverse radians !!
  v2      : 0d0  ,$   v in inverse radians !!
  baseline2: 0d0 ,$   Projected Baseline Length in Meters
  pa2      : 0d0 ,$   Position Angle of Baseline (Degrees E of N)
  u3coord : 0d0 ,$	; Closing UV coordinate (derived from uv1,uv2)
  v3coord : 0d0 ,$
  u3      : 0d0  ,$   u in inverse radians !!
  v3      : 0d0  ,$   v in inverse radians !!
  baseline3: 0d0 ,$   Projected Baseline Length in Meters
  pa3      : 0d0 ,$   Position Angle of Baseline (Degrees E of N)
  sta_index: [1,2,3],$
  flag   : byte(0) $
  }

end
