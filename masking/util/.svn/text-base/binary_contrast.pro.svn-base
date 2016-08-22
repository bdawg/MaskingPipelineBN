;;This procedure does least-squares fitting to closure-phase and
;;visibilty amplitude over a grid. It doesn't seem to work for
;;visibility amplitude.

function binary_contrast,  infile,  printfile = printfile,  init_crat = init_crat,  usev2 = usev2

common t3block,  t3data,  apriori, vis2data, usevis,  cp_cinv,  proj, proj_err
usevis = 0
cp_cinv = -1
proj = -1
p = ''

if not keyword_set(init_crat) then init_crat = 1e4

;;Read-in t3data and v3data
extract_t3data,  file = infile,  t3data
extract_vis2data,  file = infile,  vis2data
;;t3data.t3phi /=  3 ;Eliminate extra errors as a test...

;;!!! Multiply errors by 3, so that they are realistic in the context
;;of the average of 9 calibrator files.
vis2data.vis2err *= 3 
t3data.t3phierr *= 3

nbsdf = n_elements(vis2data) -(n_elements(t3data)/float(n_elements(vis2data))*3.+1) 
nbs = n_elements(t3data)
nv2 = n_elements(vis2data) 

;;Find the extra errors to add to Closure-phase to get a chi^2 of 1.0.
extra_error =  1e-3 * 10.^(alog10(1e4)*findgen(40)/40.)
newchi2 =  dblarr(40)
for i =  0, 39 do newchi2[i] =  total(t3data.t3phi^2/(t3data.t3phierr^2+extra_error[i]^2)) 
newchi2 /= nbs 
if (newchi2[0] gt 1) then begin
 extra_error = interpol(extra_error,  newchi2,  1)
 print,  extra_error,  format = '("Extra Closure-phase error: ",F6.3)'  
endif else begin
 print,  'Chi2 for closure-phase: ',  newchi2[0]
 extra_error = 0
endelse
;;Add the extra error, then multiply by an appropriate constant to
;;take into account degeneracies between closure triangles.
t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)*sqrt(nbs/nbsdf) ;+ 10.0

;;Find the extra errors to add to V^2 to get a chi^2 of 1.0.
extra_error =  1e-5 * 10.^(alog10(1e4)*findgen(40)/40.)
newchi2 =  dblarr(40)
;expr =  'p[0] + p[1]*X'
;p =  mpfitexpr(expr, r,  vis2data.vis2data/visib^2,  vis2data.vis2err/visib^2, [1.0, 0],  /quiet,  yfit = atm)
for i =  0, 39 do newchi2[i] =  total((vis2data.vis2data-1)^2/(vis2data.vis2err^2+extra_error[i]^2)) 
newchi2 /= nv2 
if (newchi2[0] gt 1) then begin
 extra_error = interpol(extra_error,  newchi2,  1)
 print,  extra_error,  format = '("Extra V^2 error: ",F8.5)'  
endif else begin
 print,  'Chi2 for V^2: ',  newchi2[0]
 extra_error = 0.0
endelse
if not keyword_set(usev2) then begin
 print,  'Adding 0.2 to V^2 errors (i.e. not using V^2)'
 extra_error += 0.2
endif else usevis = 1
vis2data.vis2err = sqrt(vis2data.vis2err^2 + extra_error^2)

;;Set up a grid of parameters
r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./8/maxr)
;;The next lines are a little hacked-up, but works for the 9h mask,
;;where the baseline separation really defines the field-of-view...
s =  sort(r) 
if (keyword_set(maxsep) eq 0) then maxsep =  rad2mas(1./2./min([r, sqrt(!pi*r[s[5]]^2/12.)]))
print,                     minsep, maxsep, $
  format = '("Min and max separations in search (mas): ", 2F7.1)'
if (p ne '') then printf,1,minsep, maxsep, $
  format = '("Min and max separations in search (mas): ", 2F7.1)'
;Set the minimum angle to be 1 radian of fringe phase at the max
;baseline...
delsep =  rad2mas(0.7/2./!pi/maxr)*1.5 ;0.7 radians. This is kind-of arbitrarily chosen...
nr   =  round((maxsep-minsep)/delsep) ;+2
delang = (delsep/maxsep)*180./!pi
nang =  round( 360./delang )
delang = 360./nang 
params =  dblarr(5, nr, nang)
params[3, *, *] = 0.1
params[4, *, *] = 0.1
for i = 0, nr-1 do for j = 0, nang-1 do begin
; if (i le 3) then  params[0, i, j] = minsep+i*delsep/2. $;separation: finer sampling
; else params[0, i, j] = minsep+(i-2)*delsep ;separation
 params[0, i, j] = minsep+i*delsep ;separation
 params[1, i, j] = j*delang        ;position angle
 params[2, i, j] = init_crat       ;brightness ratio
endfor
apriori = {val:double([0, 0, 100]),err:double([1e3, 1e3, 1e7])}

;;Now go through a grid and find the best solutions...
contrast_arr =  dblarr(nr, nang)
sig_arr = contrast_arr
for i =  0, nr-1 do for j = 0, nang-1 do begin
  apriori.val[0:1] = params[0:1, i, j]
  p1 =  params[*, i, j]
  p0 = p1
  p2 = p1
  modelt3 = binary_t3data(p1,t3data=t3data)
 ; modelv2 = binary_vis2atm(p1,vis2data=vis2data,  atm = atm,  /assym)
  modelv2 = binary_vis2data(p1,vis2data=vis2data)
  ;;This is a simple linear approximation to the contrast.
  contrast_arr[i, j] = 1./params[2, i, j]*$
   total([t3data.t3phi*modelt3.t3phi/t3data.t3phierr^2, (1.0-vis2data.vis2data)*(1.0-modelv2.vis2data)/vis2data.vis2err^2])/$
   total([modelt3.t3phi^2/t3data.t3phierr^2, (1.0-modelv2.vis2data)^2/vis2data.vis2err^2])
  ;;Now we can do a pretty dumb thing. See what the chi^2 would be at
  ;;a contrast of +/- 1e-6 from this ratio. This clearly _doesn't_
  ;;                  quite work for vis2atm
;  chi2s =  dblarr(21)
;  for k = 0, 20 do begin
;    p0[2] = 1/(contrast_arr[i, j]-1e-6*(k -10))
;    modelt3 = binary_t3data(p0, t3data = t3data)
;  modelv2 = binary_vis2atm(p0,vis2data=vis2data,  atm = atm)
;;    modelv2 = binary_vis2data(p0, vis2data = vis2data)
;    chi2s[k] =  total([(t3data.t3phi-modelt3.t3phi)^2/t3data.t3phierr^2, $
;                       (modelv2.vis2data-vis2data.vis2data)^2/vis2data.vis2err^2],  /double)
;  endfor
;  if (i eq 2) then begin
;    plot,  chi2s-min(chi2s)    
;    stop
;  endif
  bestp = params[*, i, j]
  bestp[2] = 1./contrast_arr[i, j]
  apriori.err[0:1] = 1.0 ;;An attempt to get this to work with better errors. Doesn't seem to do anything.
  dummy =  binary_oifits_chi2(bestp, sig = sig,  corr = cor,  covar = covar)
  apriori.err[0:1] = 1e3
  sig_arr[i, j] =  sig[2]
endfor
sig_arr2 =  sig_arr*contrast_arr^2
mnsig = total(sig_arr2,2)/(size(sig_arr2))[2]
mncontrast =  total(contrast_arr,2)/(size(contrast_arr))[2]
sig5_old =  mncontrast + 5*mnsig
sig5 = contrast_arr + 5*sig_arr2
sig5 = total(sig5, 2)/(size(sig5))[2]
sep =  reform(params[0, *, 0])
plot, sep, -2.5*alog10(sig5)
print,  sep,  format = '('+strtrim(nr, 2)+'F6.1 )'
print,  -2.5*alog10(sig5),  format = '('+strtrim(nr, 2)+'F6.1)'

return,  sig5

end

