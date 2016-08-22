; This program generates index arrays for an N-hole mask.
; 
; INPUT: 	n_holes 	- number of holes in the array
; OUTPUT	n_baselines
;		n_bispect
;		n_cov
;		h2bl_ix		- holes to baselines index
;		bl2h_ix		- baselines to holes index
;		bs2bl_ix	- bispectrum to baselines index
;		bl2bs_ix	- baselines to bispectrum index
;		bscov2bs_ix	- bispectrum covariance to bispectrum index
; 
; Written (from calc_bispec_template MJI)			PGT  10Dec03 

pro index_mask,n_holes,n_baselines,n_bispect,n_cov, $
    h2bl_ix,bl2h_ix,bs2bl_ix,bl2bs_ix,bscov2bs_ix, nobscov=nobscov

n_holes = long(n_holes)
n_baselines = n_holes*(n_holes-1)/2
n_bispect = n_holes*(n_holes-1)*(n_holes-2)/6
n_cov = n_holes*(n_holes-1)*(n_holes-2)*(n_holes-3)/4
;should be : ?? n_cov = n_holes*(n_holes-1)*((n_holes-2)^2)/4

; Make ordering vectors:

; Given a pair of holes i,j h2bl_ix(i,j) gives the number of the baseline
h2bl_ix = intarr(n_holes,n_holes)
count = 0
for i = 0,n_holes-2 do for j = i+1,n_holes-1 do begin
 h2bl_ix[i,j] = count
 count = count+1
endfor

; Given a baseline, bl2h_ix gives the 2 holes that go to make it up
bl2h_ix = intarr(2,n_baselines)
count = 0
for i = 0,n_holes-2 do for j = i+1,n_holes-1 do begin
 bl2h_ix[0,count] = i
 bl2h_ix[1,count] = j
 count = count +1
endfor

; Given a point in the bispectrum, bs2bl_ix gives the 3 baselines
;    which make the triangle.
; bl2bs_ix gives the index of all points in the bispectrum 
;    containing a given baseline
bs2bl_ix = intarr(3,n_bispect)
temp = intarr(n_baselines) ;N_baselines * a count variable
print, 'Indexing bispectrum...'
bl2bs_ix = intarr(n_baselines,n_holes-2)
count = 0
for i = 0,n_holes-3 do for j = i+1,n_holes-2 do for k = j+1,n_holes-1 do begin
 bs2bl_ix[0,count] = h2bl_ix[i,j]
 bs2bl_ix[1,count] = h2bl_ix[j,k]
 bs2bl_ix[2,count] = h2bl_ix[i,k]
 bl2bs_ix[bs2bl_ix[0,count],temp[bs2bl_ix[0,count]]] =count
 bl2bs_ix[bs2bl_ix[1,count],temp[bs2bl_ix[1,count]]] =count
 bl2bs_ix[bs2bl_ix[2,count],temp[bs2bl_ix[2,count]]] =count
 temp[bs2bl_ix[0,count]] = temp[bs2bl_ix[0,count]]+1
 temp[bs2bl_ix[1,count]] = temp[bs2bl_ix[1,count]]+1
 temp[bs2bl_ix[2,count]] = temp[bs2bl_ix[2,count]]+1 
 count = count+1
endfor

;Now for the long part - index the bispectral covariance...
if (keyword_set(nobscov) eq 0) then begin 
 print, 'Indexing the bispectral covariance...'
 bscov2bs_ix = intarr(2,n_cov)
 count = long(0)
 for i = 0,n_bispect-2 do for j = i+1,n_bispect-1 do begin
  if ((bs2bl_ix[0,i] eq bs2bl_ix[0,j]) or (bs2bl_ix[1,i] eq bs2bl_ix[0,j]) or (bs2bl_ix[2,i] eq bs2bl_ix[0,j]) or $
      (bs2bl_ix[0,i] eq bs2bl_ix[1,j]) or (bs2bl_ix[1,i] eq bs2bl_ix[1,j]) or (bs2bl_ix[2,i] eq bs2bl_ix[1,j]) or $
      (bs2bl_ix[0,i] eq bs2bl_ix[2,j]) or (bs2bl_ix[1,i] eq bs2bl_ix[2,j]) or (bs2bl_ix[2,i] eq bs2bl_ix[2,j]) ) then begin
   bscov2bs_ix[0,count] = i
   bscov2bs_ix[1,count] = j
   count = count + 1;
  endif
 endfor
endif

end
