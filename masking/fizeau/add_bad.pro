;######## RELEASE CLIP HERE ######## 
;This procedure adds bad baselines based on bad holes, and adds bad
;bispect based on bad baselines

pro add_bad,  bad_holes,  bad_baselines,  bad_bispect,  bl2h_ix,  bs2bl_ix,  good_baselines,  good_bispect

;Create bad baselines from bad holes
if (bad_holes[0] ne -1) then for i = 0,n_elements(bad_holes)-1 do begin
  new_bad = where((bl2h_ix[0,*] eq bad_holes[i]) or $
                  (bl2h_ix[1,*] eq bad_holes[i]))
  if (bad_baselines[0] eq -1) then bad_baselines = new_bad else $
   bad_baselines = [bad_baselines,new_bad]
endfor
bad_baselines = bad_baselines[uniq(bad_baselines,sort(bad_baselines))]

;Create a bad vector for the bispectrum...
if bad_baselines[0] ne -1 then for i = 0,n_elements(bad_baselines)-1 do begin
 new_bad = where((bs2bl_ix[0,*] eq bad_baselines[i]) or $
                 (bs2bl_ix[1,*] eq bad_baselines[i]) or $
                 (bs2bl_ix[2,*] eq bad_baselines[i]))
 if (bad_bispect[0] eq -1) then bad_bispect = new_bad else $
  bad_bispect = [bad_bispect,new_bad]
endfor
bad_bispect = bad_bispect[uniq(bad_bispect,sort(bad_bispect))]

n_baselines = (size(bl2h_ix))[2]
good_baselines = indgen(n_baselines)
if (bad_baselines[0] ne -1) then begin
 good_baselines[bad_baselines] = -1
 good_baselines = good_baselines(where(good_baselines ne -1))
endif

n_bispect = (size(bs2bl_ix))[2]
good_bispect = indgen(n_bispect)
if (bad_bispect[0] ne -1) then begin
 good_bispect[bad_bispect] = -1
 good_bispect = good_bispect(where(good_bispect ne -1))
endif

end
