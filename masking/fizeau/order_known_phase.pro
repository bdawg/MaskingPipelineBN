;######## RELEASE CLIP HERE ######## 
;This script will find ordering vectors for bl2h_ix for easy
;chi^2 calculations.

pro order_known_phase, bl2h_ix, known_ph, known_dir

n_holes = max(bl2h_ix) + 1
n_baselines = (size(bl2h_ix))[2]

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
if (ntry eq 1000) then return
if (num_ordered ne n_elements(known_ph)) then goto, try_again
known_dir = known_dir[known_order]
known_ph = known_ph[known_order]

end
