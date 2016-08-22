;;Removes problems with our oifits files, to make them
;;compatible with the oifits conventions.
;;
;;This will fix:
;;          - inconsistently defined arrname and insname entries
;;          - duplicated tables
;;          - problems with veltyp and veldef in oitarget            
;;
;; Problems:
;;          - This will erase any extra tables added to the oifits file!
;;          - The "flag" fields are not fixed here, due to IDL not having
;;             a logical or boolean type (True/False). It reads "False"
;;             as 70 though, so this will set flag=70 in the hope it works both ways
;;
;;
;;Use /flag to set flag=70 in all tables (which is what False reads into IDL as)
;;Use /python to make oifits compatible with python's oifits module (breaks binary_grid)
;; Author ACC May '13
pro fix_oifits,oifitsin,oifitsout,flag=flag,python=python

if not keyword_set(flag) then flag=0 else flag=70
if not keyword_set(target) then target=''

read_oidata, oifitsin, oiarray, oitarget,oiwavelength,oivis,oivis2, oit3

;;redo oiarray:
if keyword_set(oiarray) then begin
   oiarray=oiarray[0]
   oiarray.arrname =oit3[0].arrname
endif

;;redo oitarget
oitarget=oitarget[0]
oitarget.veltyp='UNKNOWN'
oitarget.veldef='OPTICAL'

;;redo oiwavelength
oiwavelength=oiwavelength[0]

;;and target id
for j=0,n_elements(oit3.target_id)-1 do oit3[j].target_id=0

;;and t3
for j=0,n_elements(oit3)-1 do begin
   if keyword_set(oiarray) then oit3[j].arrname=oiarray[0].arrname
   oit3[j].insname=oiwavelength[0].insname
   if keyword_set(python) then oit3[j].sta_index=0 ;;to follow standard we need sta_index=0, but this breaks binary_grid.
   oit3[j].target_id=0
endfor

;;and vis2
for j=0,n_elements(oivis2)-1 do begin
    if keyword_set(oiarray) then oivis2[j].arrname=oiarray[0].arrname
   oivis2[j].insname=oiwavelength[0].insname
   oivis2[j].target_id=0
   if keyword_set(python) then oivis2[j].sta_index=0 ;;to follow standard we need sta_index=0, but this breaks binary_grid.
endfor

;;and change "flag" since everyone uses it for different purposes.
;; Note that technically we want flag=False, but IDL has no equivalent
;; data type. IDL reads False as 70 for some reason, so maybe it works
;; both ways?
if keyword_set(flag) then begin
   for j=0,n_elements(oit3.flag)-1 do *oit3[j].flag=flag
   for j=0,n_elements(oivis2.flag)-1 do *oivis2[j].flag=flag
endif

write_oidata, oifitsout, oiarray, oitarget, oiwavelength, 0, oivis2, oit3
end
