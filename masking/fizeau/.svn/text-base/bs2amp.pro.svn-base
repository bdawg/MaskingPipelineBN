;######## RELEASE CLIP HERE ######## 
;This procedure converts a bispectrum (input as a file) to complex
;visibilties. It will only work for medium-high S/N, where at least
;roughly 2*n_holes baselines have integrated V^2 S/N > 4.
;Also, it didn't always seem to converge, so I've intuitively
;changed the fixed-point algorithm, but this might be dodgy.
;Keywords:
;diagonal_cov:  Set this to assume a diagonal covariance matrix.

pro bs2amp, mf_filestring, v2, v2_cov, bs, bs_var, bs_cov, diagonal_cov=diagonal_cov, amplitude=amplitude, amp_err=amp_err, quiet=quiet

hsn = 8.0 ;High signal-to-noise constant
;Restore the file containing the bispectrum...
restore, mf_filestring
;Create a (large!) Bispectral amplitude covariance matrix
bsacov =fltarr(n_bispect,n_bispect) 
bsacov[bscov2bs_ix[0,*],bscov2bs_ix[1,*]] = bs_cov[0,*]
bsacov = bsacov + rotate(bsacov,4)
x = indgen(n_bispect)
bsacov[x,x] = bs_var[0,*]

;Initialise some variables
x = sqrt(u^2+v^2)
in = sort(x)

;Step 1: Find initial estimates for amplitude based on V^2
dummy = cov2cor(v2_cov,sig=v2sig)
t = ((v2/v2sig - 1.0) < 1.0)> 0.0 ;Call the intermediate S/N regime between 1 and 2 sigma from 0
amp = t*sqrt(v2>0)
amp_err = t*amp*0.5*v2sig/v2 + $ ;obvious high S/N formula
  (1.0-t)*sqrt((v2>0)+2*v2sig)/2.0 ;;match the 2-sigma max V^2 for low S/N


;Step 2: Use the covariance matrix for bispectral amplitude to find
;better estimates for amplitude...
; Note that the covariance between V^2 and bispectral amplitude is not
; used here - the V^2 estimate for amplitude and the mean bispectral
; estimate are implicitly assumed independent. A formula for this in
; the case of only one bispectral estimate is:
; Cov(sqrt(v2), bs/a1/a2) = Cov(v2, bs)/2/a1/a2/a0, with a0 our best
; estimate for the amplitude of the current baseline. est and the
; est_cov matrix could be enlarged appropriately.
trynum = 0
sn = amp/amp_err
new_hsn_num = n_elements(where(sn gt hsn))
if not keyword_set(quiet) then print, new_hsn_num/float(n_baselines)*100.0, ' per cent of baselines have high S/N'
try_again:
hsn_num = new_hsn_num
trynum = trynum + 1
if not keyword_set(quiet) then print, 'Trying to find amplitude, iteration ',strtrim(trynum,2)
for tempi = 0,n_baselines-1 do begin
 i = in[tempi]
;Find bispectral points where the two other baselines have S/N greater
;than hsn[=8.0] (otherwise errors are seriously skewed in sqrt and divisions)
  w = where(((bs2bl_ix[0,*] eq i) and (sn[bs2bl_ix[1,*]] gt hsn) and (sn[bs2bl_ix[2,*]] gt hsn)) or $
            ((bs2bl_ix[1,*] eq i) and (sn[bs2bl_ix[0,*]] gt hsn) and (sn[bs2bl_ix[2,*]] gt hsn)) or $
            ((bs2bl_ix[2,*] eq i) and (sn[bs2bl_ix[0,*]] gt hsn) and (sn[bs2bl_ix[1,*]] gt hsn)))
  if (w[0] ne -1) then begin ;Now we can find a bispectral estimate for this new baseline's amplitude
   n_est = n_elements(w)
   a1 = fltarr(n_est)
   a2 = fltarr(n_est)
   a1err = fltarr(n_est)
   a2err = fltarr(n_est)
   one_vect =  replicate(1.0, n_est)
   for j = 0,n_est-1 do begin
    if (bs2bl_ix[0,w[j]] eq i) then begin
     a1[j] = amp[bs2bl_ix[1,w[j]]]
     a1err[j] = amp_err[bs2bl_ix[1,w[j]]]
     a2[j] = amp[bs2bl_ix[2,w[j]]]
     a2err[j] = amp_err[bs2bl_ix[2,w[j]]]
    endif
    if (bs2bl_ix[1,w[j]] eq i) then begin
     a1[j] = amp[bs2bl_ix[0,w[j]]]
     a1err[j] = amp_err[bs2bl_ix[0,w[j]]]
     a2[j] = amp[bs2bl_ix[2,w[j]]]
     a2err[j] = amp_err[bs2bl_ix[2,w[j]]]
    endif
    if (bs2bl_ix[2,w[j]] eq i) then begin
     a1[j] = amp[bs2bl_ix[1,w[j]]]
     a1err[j] = amp_err[bs2bl_ix[1,w[j]]]
     a2[j] = amp[bs2bl_ix[0,w[j]]]
     a2err[j] = amp_err[bs2bl_ix[0,w[j]]]
    endif
   endfor
   ;Find the estimators and their covariance matrix
   est = abs(bs[w])/a1/a2
   est_cov = fltarr(n_est,n_est)
   for j = 0,n_est-1 do for k = 0,n_est-1 do $
    est_cov[j,k] = bsacov[w[j],w[k]]/a1[j]/a1[k]/a2[j]/a2[k]
   for j = 0,n_est-1 do $
    est_cov[j,j] = est_cov[j,j] + est[j]^2*( (a1err[j]/a1[j])^2 + (a2err[j]/a2[j])^2)
   ;Use the Graybill-Deal composite estimator
   if (n_est gt 1) then begin
    invcov =  invert(est_cov)
    kvect = transpose(one_vect)#invcov/(transpose(one_vect)#invcov#one_vect)[0]
    newamp =  kvect#est
    sdev =  sqrt(kvect#est_cov#transpose(kvect))
   endif else begin
     sdev = sqrt(est_cov)
     newamp = est
   endelse
   if (sn[i] gt 0.5*hsn) then begin ;Include the straight V^2 estimator
    v2amp = sqrt(v2[i])
    v2amperr = v2amp*0.5*v2sig[i]/v2[i]
    amp[i] =  wtmn([newamp,v2amp], [sdev, v2amperr], new_sdev) 
    amp_err[i] = new_sdev
   endif else begin ;In this case only the bispectral estimator is good enough
    amp[i] = newamp
    amp_err[i] = sdev
   endelse
 endif
endfor
new_hsn_num = n_elements(where(amp/amp_err gt hsn))
if not keyword_set(quiet) then print, new_hsn_num/float(n_baselines)*100.0, ' per cent of baselines have high S/N'
;plothist, amp/amp_err
;wait, 0.5 ;Let me see this histogram...
;I'm not quite sure how many times this loop should be run, but it
;should certainly go until all possible baselines have 'high' S/N
;And no way it should be run more than ~20 times
if (new_hsn_num ne hsn_num and trynum lt 20) then goto, try_again

amplitude = amp

end
