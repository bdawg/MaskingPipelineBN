;######## RELEASE CLIP HERE ######## 
 ; This function defines the olog data structure for use
; in calibrating data. It is filled with appropriate
;numbers either manually or in an inquire_ program.
;Setting the keyword clog will reset bad_holes etc. variables, unless
;new ones are specified.

function make_clog,  olog,  clog = oldlog,  bad_holes = bad_holes,  bad_baselines = bad_baselines,  bad_bispect = bad_bispect

if not keyword_set(bad_holes) then bad_holes = -1
if not keyword_set(bad_baselines) then bad_baselines = -1
if not keyword_set(bad_bispect) then bad_bispect = -1
clog = {correction_const:0.0, nows:0, apply_phscorr:1, subtract_cp:1, diagonal_cov:1,  $
  bad_holes:bad_holes,  bad_baselines:bad_baselines, bad_bispect:bad_bispect,  cal4src:olog.cal4src,  average_src:0,  $
primary_oifits:'',  outname:'',  v2div:0.0}
if keyword_set(oldlog) then clog = {correction_const:oldlog.correction_const, nows:oldlog.correction_const, $
  apply_phscorr:oldlog.apply_phscorr, subtract_cp:oldlog.subtract_cp, diagonal_cov:oldlog.diagonal_cov,  $
  bad_holes:bad_holes,  bad_baselines:bad_baselines, bad_bispect:bad_bispect,  cal4src:oldlog.cal4src,  average_src:oldlog.average_src,  $
primary_oifits:oldlog.primary_oifits,  outname:oldlog.outname, v2div:oldlog.v2div}

return,  clog

end
