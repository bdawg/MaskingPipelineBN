;######## RELEASE CLIP HERE ######## 
;This function fills in the olog.source_name argument based on a
; Keck Catalog.

pro fill_source_name, olog,  starlist = starlist,  max_dist = max_dist

for i = 0, n_elements(olog.source_name)-1 do $
 olog.source_name[i] = find_tname(olog.ra[i],  olog.dec[i],  $
  starlist = starlist,  max_dist = max_dist)

end
