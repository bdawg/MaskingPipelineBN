;This function returns the residuals of phase from a non-redundant mask,
;where the 'hole' phases are allowed to be modified. This can be used
;by mpfit.pro
;
; mph:   Model phases
; ph:    Measured phases
; pherr: Errors on measured phases

function nr_phase_resid, mph, ph, pherr, bl2h_ix

n_holes = max(bl2h_ix) + 1
n_baselines = (size(bl2h_ix))[2]
mvis = exp(complex(0,mph))
vis = exp(complex(0,ph))
; hole_phasor[b_order_2[known_phases[i]]] = hole_phasor[b_order_1[known_phases[i]]]
;             *model_vis[known_phases[i]]/abs(model_vis[known_phases[i]])
;This needs 1 sqrt per hole, and * is complex multiplication. Note that this assumes that the b_order
;vectors will work with the known phases (eg both in ascending order of hole_num). Once this is done,
;we can form new model visibilities thus:
; new_model_vis = model_vis*hole_phasor[b_order_2]/hole_phasor[b_order_1]
;               = model_vis*hole_phasor[b_order_2]*conj(hole_phasor[b_order_1])
known_ph = where(pherr eq 0, complement=unknown)
known_dir = intarr(n_elements(known_ph))
known_order = intarr(n_elements(known_ph))

;-----------------------------------------------------------------
;Now the tricky part, we have to join the known baselines together to
;form a web from an arbitrary starting point
holes_reached = intarr(n_holes)
holes_reached[0] = 1 
num_ordered = 0
ntry = 0
try_again:

h = (where(holes_reached eq 1))
for i = 0,n_elements(h) - 1 do begin
                                ;First find the known baselines that
                                ;go outwards from this known hole...
 b = where(bl2h_ix[0,known_ph] eq h[i])
 ;Only count those where the other hole making up the baseline is unknown...
 if (b[0] ne -1) then w = where(holes_reached[bl2h_ix[1,known_ph[b]]] eq 0) else w = -1
 if (w[0] ne -1) then begin
  b = b[w]
  if (b[0] ne -1) then for j = 0,n_elements(b)-1 do begin
    known_dir[b[j]] = 0
    known_order[num_ordered] = b[j]
    num_ordered = num_ordered+1
    holes_reached[bl2h_ix[1,known_ph[b[j]]]] = 1  
  endfor
 endif
 ;Now find the known baselines that go inwards from this known hole...
 b = where(bl2h_ix[1,known_ph] eq h[i])
 ;Only count those where the other hole making up the baseline is unknown...
 if (b[0] ne -1) then w = where(holes_reached[bl2h_ix[0,known_ph[b]]] eq 0) else w = -1
 if (w[0] ne -1) then begin
  b = b[w]
  if (b[0] ne -1) then for j = 0,n_elements(b)-1 do begin
   known_dir[b[j]] = 1
   known_order[num_ordered] = b[j]
   num_ordered = num_ordered+1
   holes_reached[bl2h_ix[0,known_ph[b[j]]]] = 1  
  endfor
 endif
endfor

ntry = ntry + 1
if (ntry eq 1000) then stop
if (num_ordered ne n_elements(known_ph)) then goto, try_again
known_dir = known_dir[known_order]
known_ph = known_ph[known_order]

;--------------------------------------------------------
;Use this knowledge to construct hole phasors, and then to correct the
;model visibilities. All model visibilities should have zero phase at
;the known_baselines...
hole_phasors = complexarr(n_holes)
hole_phasors[*] = 1.0
for i = 0,(n_elements(known_ph) - 1) do begin
 if (known_dir[i] eq 0) then $;forwards
   hole_phasors[bl2h_ix[1,known_ph[i]]] = hole_phasors[bl2h_ix[0,known_ph[i]]]*mvis[known_ph[i]] $
 else $ ;backwards
   hole_phasors[bl2h_ix[0,known_ph[i]]] = hole_phasors[bl2h_ix[1,known_ph[i]]]*conj(mvis[known_ph[i]])
endfor
new_mvis = mvis*hole_phasors[bl2h_ix[0,*]]*conj(hole_phasors[bl2h_ix[1,*]])

;-----------------------------------------------------------
;Now find the chi^2
viserr = pherr
viserr[known_ph] = 1.0 ;This number shouldn't matter apart from errors in floating arithmetic...

return, atan(new_mvis[unknown]/vis[unknown], /phase)/viserr[unknown]

end
