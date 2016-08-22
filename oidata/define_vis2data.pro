pro define_vis2data,vis2data_unit
;+
; NAME:
;     define_vis2data
; PURPOSE:
;     Define the the VIS2DATA structure which will contain the data in 
;     the OIDATA format in a easily accessible format for by IDL programs.
;
; CALLING SEQUENCE:
;     define_vis2data, vis2data_unit
;
; INPUTS:
;     None
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;      VIS2DATA_UNIT  
;		 Will contain measurements on the interferometer complex
;		 visibility. This structure can be read by WRITE_OIDATA
;		 to write standard OI-DATA FITS files.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      Used by extract_vis2data.
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
;     v0.0 2003Jul04    J.D. Monnier    Initiated
;     v0.1 2005May30    JDM		targets
;     v0.2 2006Oct22	JDM		sfu
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
vis2data_unit = { $
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
  vis2data: 0d0 ,$ 
  vis2err : 0d0 ,$
  ucoord : 0d0	,$  
  vcoord : 0d0	,$
  u      : 0d0  ,$   u in inverse radians !!
  v      : 0d0  ,$   v in inverse radians !!
  sfu    : 0d0  ,$   Baseline length in inverse Radians.
  baseline: 0d0 ,$   Projected Baseline Length in Meters
  pa      : 0d0 ,$   Position Angle of Baseline (Degrees E of N)
  sta_index: [0,0], $ Station Index
  flag   : byte(0) $
  }

end
