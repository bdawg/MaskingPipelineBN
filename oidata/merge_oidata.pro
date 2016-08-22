pro merge_oidata, outfile=outfile, infiles=infiles,  $
  oiarray_out,oitarget_out,oiwavelength_out,oivis_out,oivis2_out,oit3_out
;+
; NAME:
;     merge_oidata
; PURPOSE:
;     To merge multiple OIDATA.FITS files into a single file with multiple
;     instances of relevant OI-tables.
;
; CALLING SEQUENCE:
;     merge_oidata, outfile='hybrid.fits', infiles=['file1.fits','file2.fits',file3.fits']
; 	or
;     merge_oidata, outfile='hybrid.fits', infiles=['file1.fits','file2.fits',file3.fits'], $
;	  oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
;
;
; INPUTS:
;     outfile - Name of file to write with the merged data. If outfile is not specified
;  	        then no file is written.
;     infiles    - string array containing names of (OIDATA) fits files to merge.
;
;  Note: By default this program will keep all the information from the
;	 individual files in separate OI-TABLES so absolutely no information 
; 	 is lost (and the process can be reversed). This may necessitate
;        changing the 'INSNAME' and 'ARRNAME' to make sure they are unique.
;	 This will done automatically by the program adding numbers to common names.
;	 Any name changes will be reported.
;
; OPTIONAL INPUTS:
;     None 
;
; KEYWORDS:
;     None
;
; OUTPUTS:
;     See header of 'read_oidata' for details on each of the OIDATA Tables.
;
;      OIARRAY   Contains information on the interferometer array geometry.
;      OITARGET  Contains information on the observed targets.
;   OIWAVELENGTH Contains information on the observing bandpasses.
;      OIVIS     Contains complex visibility measurements.
;      OIVIS2    Contains visibility-squared measurements.
;      OIT3      Contains measurements of the bispectrum (triple amp and phase).
;
; 
; RESTRICTIONS:
;      None.
;
; PROCEDURE:
;      Simple issue the command, as shown in the calling sequence above.
;      This procedure will read in multiple files and write a big oidata output file.
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
;	See 'Calling Sequence' above for an example of use of this command.
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
;     v0.0 2003Jul01    J.D. Monnier    Written for OI-DATA Revision 1
;     v0.1 2006Mar03	J.D. Monnier	Bugfix in merging OI_VIS tables
;
;     To do:   Better Error checking.
;	       Way to join together data more cleverly (to share arrays or
;		  wavelength tables, etc). E.g. interactive way to
;	  	  identify targets entries which are identical and to
;		  create a virtual INSNAME detector which combine capabilities
;		  of multiple arrays, etc.
;
;  *** Please write monnier@umich.edu with bug reports. ***
;-
;---------------------------------------------------------------------------


num_infiles = n_elements(infiles)

; Read the first one
read_oidata, infiles(0), oiarray,oitarget,oiwavelength,oivis,oivis2,oit3

for i=1,num_infiles-1 do begin
  read_oidata, infiles(i), oiarray0,oitarget0,oiwavelength0,oivis0,oivis20,oit30

;Get lists of targets, insname, arrays

; First combine OI Target. Easiest since only one table allowed (and required) per file.

; Before concatenating, one must check uniqueness issues!
; do targets first since target arrays must exist in each file.

if (n_elements(oitarget0) eq 0 or n_elements(oitarget) eq 0) then begin
   print, 'There is supposed to be a OITARGET table in each file'
   stop
endif

for j=0,n_elements(oitarget0)-1 do begin
  ; only target_IDs must be uniq.. not Target names.... Assume there are no
  ; inconsistencies in the files that are being merged (assume no target_ids repeated)
    repeat_in = where(oitarget.target_id eq  oitarget0(j).target_id,ct_repeat)
    if (ct_repeat ge 1) then begin
      ; Find lowest available target_id in both oitarget,oitarget0
      candidate_id=0
      idlist = [oitarget.target_id, oitarget0.target_id]
      repeat_id = where( idlist eq candidate_id,ct_id )
      while (ct_id ge 1) do begin
        candidate_id=candidate_id +1
        repeat_id = where( idlist eq candidate_id,ct_id)
      endwhile
  ; ok.. now we have a new candidate_ID -- so we must replace the old ID with the new ID
  ; in the target table and all data tables (that exist).
    old_id = oitarget0(j).target_id
    oitarget0(j).target_id = candidate_id
    if (n_elements(oivis0) ge 1) then begin
      in0=where(oivis0.target_id eq old_id,ct)
      if (ct gt 0) then oivis0(in0).target_id = candidate_id
    endif
   if (n_elements(oivis20) ge 1) then begin
      in0=where(oivis20.target_id eq old_id,ct)
      if (ct gt 0) then oivis20(in0).target_id = candidate_id
    endif 
   if (n_elements(oit30) ge 1) then begin
      in0=where(oit30.target_id eq old_id,ct)
      if (ct gt 0) then oit30(in0).target_id = candidate_id
    endif
    endif ; if there is repeat
endfor

if ( n_elements(oiarray) ge 1 and n_elements(oiarray0) ge 1 ) then begin 
  
for j=0,n_elements(oiarray0)-1 do begin
  arr=strtrim(oiarray.arrname,2)
  arr0=strtrim(oiarray0.arrname,2)
  arrlist=[arr,arr0]
  repeat_in = where(arr eq arr0(j),ct_repeat)
  if (ct_repeat ge 1) then begin
   increment=0
   while (ct_repeat ge 1) do begin
     increment=increment+1
     strincrement=strtrim(string(increment),2)
     ;Create a unique name for the array of form "name_v1"
     newname = arr0(j)+"_v"+strincrement
     repeat_again = where(arrlist eq newname,ct_repeat) 
   endwhile
  ; use realname
   oldname = arr0(j)
   multiplev=where( strtrim(oiarray0.arrname,2) eq oldname)
   oiarray0(multiplev).arrname=newname
  
  if (n_elements(oivis0) ge 1) then begin
      in0=where(strtrim(oivis0.arrname,2) eq oldname,ct)
      if (ct gt 0) then oivis0(in0).arrname=newname
    endif
   if (n_elements(oivis20) ge 1) then begin
      in0=where(strtrim(oivis20.arrname,2) eq oldname,ct)
      if (ct gt 0) then oivis20(in0).arrname=newname
    endif
 if (n_elements(oit30) ge 1) then begin
      in0=where(strtrim(oit30.arrname,2) eq oldname,ct)
      if (ct gt 0) then oit30(in0).arrname=newname
    endif
  endif; if repeat
  
endfor
endif ; if both array/ array0 exist

; Last but not least.. OIWAVELENGTH/INSNAME tables
; like oitarget, all oidata fits files should have at least one of these so we
; needn't really have to check...
if (n_elements(oitarget0) eq 0 or n_elements(oitarget) eq 0) then begin
   print, 'There is supposed to be a OIWAVELENGTH table in each file'
   stop
endif

for j=0,n_elements(oiwavelength0)-1 do begin
  wav=strtrim(oiwavelength.insname,2)
  wav0=strtrim(oiwavelength0.insname,2)
  wavlist=[wav,wav0]
  repeat_in=where(wav eq wav0(j),ct_repeat)
  if (ct_repeat ge 1) then begin
   increment=0
   while (ct_repeat ge 1) do begin
     increment=increment+1
     strincrement=strtrim(string(increment),2)
     ;Create a unique name for the insname of form "name_v1"
     newname = wav0(j)+"_v"+strincrement
     repeat_again = where(wavlist eq newname,ct_repeat) 
   endwhile
  ; use realname
   oldname = wav0(j)
   multiplev=where( strtrim(oiwavelength0.insname,2) eq oldname)
   oiwavelength0(multiplev).insname=newname
   if (n_elements(oivis0) ge 1) then begin
      in0=where(strtrim(oivis0.insname,2) eq oldname,ct)
      if (ct gt 0) then oivis0(in0).insname=newname
    endif
   if (n_elements(oivis20) ge 1) then begin
      in0=where(strtrim(oivis20.insname,2) eq oldname,ct)
      if (ct gt 0) then oivis20(in0).insname=newname
    endif
 if (n_elements(oit30) ge 1) then begin
      in0=where(strtrim(oit30.insname,2) eq oldname,ct)
      if (ct gt 0) then oit30(in0).insname=newname
    endif
  endif; if repeat
  
endfor

; ok. now concat the tables/structures.
;OITARGET
oitarget=concat_oitable(oitarget,oitarget0)
;OIWAVELENGTH
oiwavelength0.extver=oiwavelength0.extver + max(oiwavelength.extver)
oiwavelength=concat_oitable(oiwavelength,oiwavelength0)
;OIARRAY
if (n_elements(oiarray) gt 0 and n_elements(oiarray0) gt 0) then begin
  oiarray0.extver = oiarray0.extver + max(oiarray.extver)
  oiarray=concat_oitable(oiarray,oiarray0)
endif 
if (n_elements(oiarray) eq 0 and n_elements(oiarray0) gt 0) then $
   oiarray=oiarray0
;OIVIS
if (n_elements(oivis) gt 0 and n_elements(oivis0) gt 0) then begin
  oivis0.extver = oivis0.extver + max(oivis.extver)
  oivis=concat_oitable(oivis,oivis0)
endif 
if (n_elements(oivis) eq 0 and n_elements(oivis0) gt 0) then $
  oivis=oivis0
;OIVIS2
if (n_elements(oivis2) gt 0 and n_elements(oivis20) gt 0) then begin
  oivis20.extver = oivis20.extver + max(oivis2.extver)
  oivis2=concat_oitable(oivis2,oivis20)
endif 
if  (n_elements(oivis2) eq 0 and n_elements(oivis20) gt 0) then $
  oivis2=oivis20
;OIT3
if (n_elements(oit3) gt 0 and n_elements(oit30) gt 0) then begin
  oit3.extver = oit3.extver + max(oit3.extver)
  oit3=concat_oitable(oit3,oit30)
endif 
if (n_elements(oit3) eq 0 and n_elements(oit30) gt 0) then $
  oit3=oit30

endfor


if (n_elements(outfile) ne 0) then begin
  write_oidata, outfile,oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
endif

oitarget_out=oitarget
oitwavelength_out=oiwavelength
if (n_elements(oiarray) gt 0) then oiarray_out=oiarray else delvarx,oiarray_out
if (n_elements(oivis) gt 0) then oivis_out=oivis else delvarx,oivis_out
if (n_elements(oivis2) gt 0) then oivis2_out=oivis2 else delvarx,oivis2_out
if (n_elements(oit3) gt 0 ) then oit3_out=oit3 else delvarx,oit3_out

end


