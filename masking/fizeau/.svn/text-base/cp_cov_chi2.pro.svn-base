;This function returns the chi-squared of a fit to an array
;unorm and vnorm are normalized (to max baseline) u and v coordinates.

function cp_cov_chi2,  p
common cp_cov_func,  n_baselines,  n_bispect,  n_cov, bs2bl_ix, bscov2bs_ix, meas_cov, modcov,  alphamat,  bsmodsq, x0, x1, x2, y0, y1, y2, cov_sig, ncalls
ncalls++
modcov = meas_cov
modcov[*] = 0.

;print,  'Calculating chi^2.. ',systime()  

;First, find the baseline-based phase covariance
;!!! NB, there is a fixed addition in the following line of 0.06
;degrees per baseline - an assumed minimum error.
bcov = diag_matrix(p[0:n_baselines-1]^2 + 1e-6) + p[n_baselines+1]^2*alphamat
modcov = bcov[x0, y0] + bcov[x0, y1] - bcov[x0, y2] $
       + bcov[x1, y0] + bcov[x1, y1] - bcov[x1, y2] $
       - bcov[x2, y0] - bcov[x2, y1] + bcov[x2, y2]
ix = indgen(n_bispect)
modcov[ix, ix] += p[n_baselines]^2/bsmodsq

;Now, convert this to a bispectrum-based phase covariance
;for i = 0, n_bispect-1 do for j = 0, n_bispect-1 do $
;for j = 0, n_bispect-1 do $
; modcov[*, j] = bcov[bs2bl_ix[0, *], bs2bl_ix[0, j]] + bcov[bs2bl_ix[0, *], bs2bl_ix[1, j]] - bcov[bs2bl_ix[0, *], bs2bl_ix[2, j]] $
;               +bcov[bs2bl_ix[1, *], bs2bl_ix[0, j]] + bcov[bs2bl_ix[1, *], bs2bl_ix[1, j]] - bcov[bs2bl_ix[1, *], bs2bl_ix[2, j]] $
;               -bcov[bs2bl_ix[2, *], bs2bl_ix[0, j]] - bcov[bs2bl_ix[2, *], bs2bl_ix[1, j]] + bcov[bs2bl_ix[2, *], bs2bl_ix[2, j]]

;for i = 0, n_cov-1 do begin
; bs0 = bscov2bs_ix[0, i] 
; bs1 = bscov2bs_ix[1, i]
; w0 =  where(bs2bl_ix[*, bs0] eq bs2bl_ix[0, bs1] or $
;             bs2bl_ix[*, bs0] eq bs2bl_ix[1, bs1] or $
;             bs2bl_ix[*, bs0] eq bs2bl_ix[2, bs1])
; w1 =  where(bs2bl_ix[*, bs1] eq bs2bl_ix[0, bs0] or $
;             bs2bl_ix[*, bs1] eq bs2bl_ix[1, bs0] or $
;             bs2bl_ix[*, bs1] eq bs2bl_ix[2, bs0])
; if (w0[0] eq -1 or w1[0] eq -1) then stop
; if ((w0[0] eq 2 and w1[0] ne 2) or (w1[0] eq 2 and w0[0] ne 2)) then sign = -1 else sign = 1
; modcov[bs0, bs1] = sign*p[bs2bl_ix[w0, bs0]]^2
;endfor
;modcov =  modcov+transpose(modcov)
;for i = 0, n_bispect-1 do $
; modcov[i, i] += p[n_baselines]^2/bsmodsq[i] ;+ p[bs2bl_ix[0, i]]^2 + p[bs2bl_ix[1, i]]^2 + p[bs2bl_ix[2, i]]^2
; modcov[i, i] = p[n_baselines]^2/bsmodsq[i] + p[bs2bl_ix[0, i]]^2 + p[bs2bl_ix[1, i]]^2 + p[bs2bl_ix[2, i]]^2

retval = total(((modcov-meas_cov)/cov_sig)^2)
;pb = p[0:n_baselines-1]^2
;mpb = mean(pb)
;Penalize baselines that have less than half the median error [Does
;this bias anything...?]

;retval += total( ( (mpb/4. - pb) > 0)/mpb) *n_baselines*4.

return,  retval

end
