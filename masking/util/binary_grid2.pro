;;This program is designed to be a cut-down and simpler version of
;;binary_grid, focussing on diagnostics rather than fancy statistics.
;;Keywords:
;;  infile: The input oifits file. We don't input the
;;    closure-phase covariance matrix to this program - only the data
;;    that is part of the oifits format.
;;  init_crat: The contrast ratio at which the grid search is done. NB
;;    in the high contrast regime, you'll get the same best
;;    solution no matter what init_crat is chosen.
;;  boxsize: The size of the grid box in milli-arcsec
;;  like: An output likelihood map.
;;  x: An output vector giving the axes of the likelihood map in mas.
;;  minsep: The resolution of the grid.
pro binary_grid2,  infile,  init_crat=init_crat, boxsize=boxsize,like=like,x=x, extra_error=extra_error

if (keyword_set(init_crat) eq 0) then init_crat = 250.
extract_t3data,  file = infile,  t3data
extract_vis2data,  file = infile,  vis2data
;;Deal with "bad" triangles in a simplistic way.
w = where(t3data.flag)

 if (w[0] ne -1) then t3data.t3phierr = 1e3
 if keyword_set(extra_error) then t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)
 extract_vis2data,  file = infile,  vis2data
 read_oidata,infile,oiarray,oitarget
 target_name=oitarget.target
 target_name=target_name[0]
 ;Enforce a maximum SNR of 4.0 on visibility data.
 vis2data.vis2err = sqrt(vis2data.vis2err^2 + (0.2*vis2data.vis2data)^2)

;;Find the number of rotations on the sky (i.e. independent cubes) by
;;finding how many times the first baseline in the file is repeated,
;;where we define "baseline" as containing uniqu stations indices.
;;WARNING: breaks for a >100 hole mask.
ix = vis2data.sta_index[0] + 100*vis2data.sta_index[1]
nrotations = n_elements(where(ix eq ix[0]))
nv2 = n_elements(vis2data)/nrotations
nt3 = n_elements(t3data)/nrotations

;Number of independent degrees of freedom in the data. The complicated
;expression at the end is (n_holes - 1)
ndf =  nrotations*(nv2 - (3*nt3/nv2+1) )

r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./8/maxr)
nsep=boxsize/minsep

;Now search for the best chi^2 over the grid...
chi2_arr =  fltarr(nsep, nsep)
x=(findgen(nsep)-nsep/2+0.5)*minsep
make_2d, x,x,xx, yy
rho = sqrt(xx^2+yy^2)
theta = atan(-xx,yy)*180/!pi
for i =  0, nsep-1 do for j = 0, nsep-1 do begin
  params = [rho[i,j], theta[i,j], init_crat,0.1,0.1]
  modelt3 = binary_t3data(params,t3data=t3data)
  cpresid = mod360(modelt3.t3phi - t3data.t3phi)
  chi2_arr[i, j] =  total(cpresid^2/t3data.t3phierr^2) 
endfor
;;Make reduced chi^2 equal to 1 by scaling errors, and make sure that
;;our final chi^2 variable has a minimum of ndf (crudely takes into
;;account linear dependence of closure-phases)
ndf =  n_elements(vis2data) -(n_elements(t3data)/float(n_elements(vis2data))*3.+1)
chi2_arr /= min(chi2_arr)
chi2_arr *= ndf
;;The likelihood map.
like = exp(-(chi2_arr-ndf)/2)
y=x
image_cont_deluxe, -like, /noc, /asp, xv=-x, yv=y, ytit='Dec (mas)', xtit='RA (mas)'
oplot, [0], [0], psym=2 
 ;image_cont, chi2_arr, /noc
 ;stop
end

