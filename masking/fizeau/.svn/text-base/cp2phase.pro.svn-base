;######## RELEASE CLIP HERE ######## 
;This procedure converts a bispectrum (input as a file) to complex
;visibilties. It will only work for medium-high S/N, where at least
;roughly 2*n_holes baselines have integrated V^2 S/N > 5.
;Also, it didn't always seem to converge, so I've intuitively
;changed the fixed-point algorithm, but this might be dodgy.
;Keywords:
;diagonal_cov:  Set this to assume a diagonal covariance matrix.

pro cp2phase, mf_filestring, cp, cperr, cpcov, diagonal_cov=diagonal_cov, phase=phase, ph_err=ph_err

;Restore the file containing the bispectrum...
restore, mf_filestring

;Initialise some variables (like closure phase!)
x = sqrt(u^2+v^2)
in = sort(x)
known_phase = intarr(n_baselines)
phase = fltarr(n_baselines)
ph_err = fltarr(n_baselines)

;Step 1: Find the shortest nholes-1 baselines that cover all holes
;and don't form any closed loops.
known_phase[in[0]] = 1
known_holes = intarr(n_holes)
known_holes[bl2h_ix[0,in[0]]] = 1
known_holes[bl2h_ix[1,in[0]]] = 1
while (total(known_phase) lt n_holes-1) do begin
 w = where(known_holes[bl2h_ix[0,in]] + known_holes[bl2h_ix[1,in]] eq 1)
 new_bline = in[w[0]]
 known_phase[new_bline] = 1
 known_holes[bl2h_ix[0,new_bline]] = 1
 known_holes[bl2h_ix[1,new_bline]] = 1
endwhile
old_known_phase = known_phase

;Step 2: Fill known_phase, phase, and ph_error vectors from the inside-out.
trynum = 0
try_again:
trynum = trynum + 1
print, 'Trying to find phase, iteration ',strtrim(trynum,2)
second_try = -1
for tempi = 0,n_baselines-1 do begin
 i = in[tempi]
 if (known_phase[i] eq 0) then begin
 ;Find bispectral points that contain 2 baselines of known (or assumed) phase.
  w = where((bs2bl_ix[0,*] eq i or bs2bl_ix[1,*] eq i or bs2bl_ix[2,*] eq i) and $
            (known_phase[bs2bl_ix[0,*]] + known_phase[bs2bl_ix[1,*]] + known_phase[bs2bl_ix[2,*]] eq 2))
                                ;For diagnostics, keep a record of
                                ;those baselines that will require 2 goes at this...
  if (w[0] eq -1) then begin
   if (second_try[0] eq -1) then second_try = [i] else $
    second_try = [second_try,i]
  endif else begin ;Now we can find an estimate for this new baseline
   n_est = n_elements(w)
   p1in = intarr(n_est)
   p2in = intarr(n_est)
   p1sign = fltarr(n_est)
   p2sign = fltarr(n_est)
   est_sign = fltarr(n_est)
   ;Loop through the closure triangle phase estimates to
   ;work out the signs etc...
   for j = 0,n_elements(w)-1 do begin
    if (bs2bl_ix[2,w[j]] eq i) then begin
     est_sign[j] = -1.0 
     p1sign[j] = 1.0
     p2sign[j] = 1.0
     p1in[j] = bs2bl_ix[0,w[j]]
     p2in[j] = bs2bl_ix[1,w[j]]
    endif else if (bs2bl_ix[0,w[j]] eq i) then begin
     est_sign[j] = 1.0
     p1sign[j] = 1.0
     p2sign[j] = -1.0
     p1in[j] = bs2bl_ix[1,w[j]]
     p2in[j] = bs2bl_ix[2,w[j]]
    endif else if (bs2bl_ix[1,w[j]] eq i) then begin
     est_sign[j] = 1.0
     p1sign[j] = 1.0
     p2sign[j] = -1.0
     p1in[j] = bs2bl_ix[0,w[j]]
     p2in[j] = bs2bl_ix[2,w[j]]
    endif
   endfor
   est = est_sign*(cp[w]-p1sign*phase[p1in]-p2sign*phase[p2in])  
   est_err = sqrt(cperr[w]^2+ph_err[p1in]^2+ph_err[p2in]^2)
   offset = est((where(est_err eq min(est_err)))[0])
   ;if (max(abs(est/est_err)) gt 8.0) then stop
   est = rad2mpipi(est - offset)
   phase[i] = wtmn(est, est_err, sdev)+offset ;I don't know if wtmn is appropriate here...
   ph_err[i] = sdev
   known_phase[i] = 1
  endelse
 endif
endfor
if (total(known_phase) lt n_baselines) then goto, try_again

;Step 3: Now do things properly: use the covariance matrix of closure
;phase to find the phase (iteratively?) using inside-out as a
;starting point...
print, 'Original rms phase error is: ',sqrt(mean(ph_err^2))*180/!pi, ' degrees'
new_phase = phase
new_ph_err = ph_err
known_phase = old_known_phase
count = 0
old_phase = 0
step3:
count = count + 1
;old_phase is for diagnostics if something goes wrong...
old_phase = phase
old_ph_err = ph_err
;Use a weighted mean of the previous 2 estimates if things aren't
;working so well... (Remove this if it hasn't recently been used)
ph_err = 0.7*new_ph_err + 0.3*ph_err
phase = new_phase + 0.3*rad2mpipi(phase - new_phase)

;Now that we have an initial estimate for all phases, we can
;fix the size of the following arrays
n_est = n_holes-2
p1in = intarr(n_est)
p2in = intarr(n_est)
p1sign = fltarr(n_est)
p2sign = fltarr(n_est)
est_sign = fltarr(n_est)
cov = fltarr(n_est,n_est)
k = indgen(n_est) ;A convenient indexing vector
one_vect =  replicate(1.0, n_est)
;Again, go roughly inside-out...
for tempi = 0,n_baselines-1 do begin
 i = in[tempi]
 ;For the phases that are fixed at 0, this loop doesn't apply:
 if (known_phase[i] eq 0) then begin
   ;Sort out signs of the estimators...
   w = where(bs2bl_ix[0,*] eq i or bs2bl_ix[1,*] eq i or bs2bl_ix[2,*] eq i)
   if (n_elements(w) ne n_est) then stop
   for j = 0,n_elements(w)-1 do begin
    if (bs2bl_ix[2,w[j]] eq i) then begin
     est_sign[j] = -1.0 
     p1sign[j] = 1.0
     p2sign[j] = 1.0
     p1in[j] = bs2bl_ix[0,w[j]]
     p2in[j] = bs2bl_ix[1,w[j]]
    endif else if (bs2bl_ix[0,w[j]] eq i) then begin
     est_sign[j] = 1.0
     p1sign[j] = 1.0
     p2sign[j] = -1.0
     p1in[j] = bs2bl_ix[1,w[j]]
     p2in[j] = bs2bl_ix[2,w[j]]
    endif else if (bs2bl_ix[1,w[j]] eq i) then begin
     est_sign[j] = 1.0
     p1sign[j] = 1.0
     p2sign[j] = -1.0
     p1in[j] = bs2bl_ix[0,w[j]]
     p2in[j] = bs2bl_ix[2,w[j]]
    endif
   endfor
   ;Make a covariance matrix of the estimators for phase[i]:
   for j = 0,n_est-1 do cov[k,j] = cpcov[w[k],w[j]]*est_sign[k]*est_sign[j]
   est_err = sqrt(cperr[w]^2+ph_err[p1in]^2+ph_err[p2in]^2)
   cov[k,k] = est_err^2
   ;Form a simple vector of phase estimators.
   est = est_sign*(cp[w]-p1sign*phase[p1in]-p2sign*phase[p2in])
   offset = est((where(est_err eq min(est_err)))[0])
   est = rad2mpipi(est - offset)
   old_est = est
   ;Now, combine these estimators in a way that includes cpcov knowledge
   ;This involves diagonalising the covariance matrix
   ;i.e. we want to find A such newest=A#est, that E(A#est) = E(est)
   ; and E(newest#newest^T) is diagonal, which means that
   ; A#cov#A^T = D, with D diagonal.
;   orig_cov = cov
;   trired, cov, err, e
;   triql, err, e, cov
   ;Now err is the eigenvalues of orig_cov, cov is the eigenvectors
   ;norms is required as we need two constraints on a.
;   norms = total(cov,1)
;   A = transpose(cov)
;   for j = 0,n_est-1 do A[j,*]=A[j,*]/norms[j]
;   est = A # est
;   err = sqrt(err/norms^2) ;Equal to the diagonal of A#orig_cov#A^T
   ;err is now the standard deviation of independent estimators est .. 
;   est = transpose(cov) # est / norms;!!!
;   err = sqrt(err/norms^2);!!!
;   new_phase[i] = wtmn(est, err, sdev)+offset

;This new method is calles the Graybill-Deal composite estimator. See
;KellerWigton.pdf. Error calculated analytically, assuming that kvect
;is not a random vector (of course it is, because cov has errors)
;This new method is calles the Graybill-Deal composite estimator. See
;KellerWigton.pdf. Error calculated analytically, assuming that kvect
;is not a random vector (of course it is, because cov has errors)
   if (keyword_set(diagonal_cov)) then new_phase[i] = wtmn(old_est, est_err, sdev, mse=mse)+offset $
   else begin
    invcov =  invert(cov)
    kvect = transpose(one_vect)#invcov/(transpose(one_vect)#invcov#one_vect)[0]
    new_phase[i] =  kvect#old_est+offset
    sdev =  sqrt(kvect#cov#transpose(kvect))
    if (sdev lt 0.2*sqrt(1./total(1./est_err^2))) then stop
   endelse
   new_ph_err[i] = sdev > sqrt(1./total(1./est_err^2))
 endif
endfor
if (count lt 5 or (count mod 10) eq 0) then begin
 print,  'Iteration: ',  count,  ' to find phase from closure phase.'
 change = sqrt(mean(rad2mpipi(phase-new_phase)^2))*180/!pi
 print, 'rms change in phase: ', change, ' degrees'
 change = sqrt(mean(rad2mpipi(phase-new_phase)^2/(ph_err>0.001)^2))
 print, 'rms normalised change in phase: ', change
endif
if (count mod 100 eq 0) then begin
  print, 'This really isn''t converging. Try setting the diagonal_cov keyword.'
  stop
  goto, step3
  endif
if (change gt 0.1) then goto, step3

ph_err = new_ph_err
phase = new_phase

;The following plot was annoying. Now removed.
;plothist, ph_err*180/!pi, bin=0.3
print, 'Final rms phase error is: ',sqrt(mean(ph_err^2))*180/!pi, ' degrees'

end
