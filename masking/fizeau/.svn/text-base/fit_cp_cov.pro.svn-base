;This function takes a closure-phase covariance matrix and fits an
;(n_baselines+2)-parameter "filtered" covariance matrix to it,
;enabling relatively error-free inversion.
;(use la_eigenproblem to return the eigenvalues and check that they
;make sense)
;restore,  'cubeinfo5Jun07.idlvar'
;src = where(total(clog.cal4src,1) gt 0)
;restore, 'bs' + olog.cube_fname[src[0], 1] + '.idlvar'

function fit_cp_cov,  cov,  bs, mf_file,  init_p = init_p,  bestp = bestp,  alphamat = alphamat_out
;common cp_cov_func,  n_baselines,  n_bispect,  n_cov, bs2bl_ix, bscov2bs_ix, meas_cov, modcov,  alphamat,  bsmodsq
common cp_cov_func,  n_baselines,  n_bispect,  n_cov, bs2bl_ix, bscov2bs_ix, meas_cov, modcov,  alphamat,  bsmodsq, x0, x1, x2, y0, y1, y2, cov_sig, ncalls

meas_cov = cov
restore, !ROOT_DIR + 'templates/'+ mf_file
bsmodsq =  modsq(bs)

;Find the covariance matrix of baseline-phase corresponding to
;tip/tilt errors. We just find the variance of the x and y modes 
;then add them.
uvmax =  max(sqrt(u^2+v^2))
unorm = u/uvmax
vnorm = v/uvmax
bx =  unorm^3
;mx =  fltarr(n_bispect)
;for i = 0, n_bispect-1 do mx[i] = bx[bs2bl_ix[0, i]] + bx[bs2bl_ix[1, i]] - bx[bs2bl_ix[2, i]]
by =  vnorm^3
;my =  fltarr(n_bispect)
;for i = 0, n_bispect-1 do my[i] = by[bs2bl_ix[0, i]] + by[bs2bl_ix[1, i]] - by[bs2bl_ix[2, i]]
;alphamat =  mx#transpose(mx) + my#transpose(my)
alphamat =  bx#transpose(bx) + by#transpose(by)

make_2d,  bs2bl_ix[0, *], bs2bl_ix[0, *],  x0, y0
make_2d,  bs2bl_ix[1, *], bs2bl_ix[1, *],  x1, y1
make_2d,  bs2bl_ix[2, *], bs2bl_ix[2, *],  x2, y2

xi = 0.1*identity(n_baselines+2)
ix =  indgen(n_bispect)
cperr = sqrt(meas_cov[ix, ix])
if keyword_set(init_p) then bestp = init_p else $
 bestp = [replicate(median(cperr), n_baselines)/2, sqrt(median(bsmodsq))/10.,  0.003]
cor =  cov2cor(cov,  sig = sig)
cov_sig = sig#transpose(sig)
cov_sig[ix, ix] /= sqrt(n_bispect) ;Make sure that the covariances are fit first.
;Amoeba: 12:00:40 start, 1000 evaluations should be 12:02:20
;Powell needs 22353 evaluations to converge at 1e-3
;Powell needs  = - 22353 evaluations to converge at 1e-2
;      284952.   13715
;      284503.   22353
ncalls=0
powell,  bestp,  xi,  3e-3,  bestchi2,  'cp_cov_chi2'
;bestp2=bestp
;bestchi2_2=bestchi2
;print, bestchi2, ncalls
;powell,  bestp,  xi,  1e-3,  bestchi2,  'cp_cov_chi2'
;print, bestchi2, ncalls
;plot, (bestp-bestp2)/median(abs(bestp))
;stop
print,  'CP Cov Chi^2: ',  bestchi2
;Now, make sure that there is enough of a random error so that the
;returned matrix can be inverted.
bestp =  abs(bestp)
bestp[n_baselines] = bestp[n_baselines] > sqrt(median(bestp^2)*median(bsmodsq)/100.)
dummy =  cp_cov_chi2(bestp)
alphamat_out = alphamat ;For diagnostics only. 1e-2 gives 8:39:15 -> more than 8:54. Whoa.

return,  modcov

end
