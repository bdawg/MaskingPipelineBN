;This fuction returns a mean and covariance matrix of a bunch of
;parameters, given some input data. Two variables, data and errors can
;be passed through this routine to the chi^2. It is of course
;important that 'data' and 'errors' are both passed by reference, ie
;without the use of subscripts...
;
;my_funct:     A string containing the name of the chi^2 calculation
;              function.
;params:       Initial parameters - must be reasonable within errors.
;step_size:    Step size for each parameter
;num_iter:     Number of iterations. I think to get a 10% error on
;              covariances you need 400*(sig/step_size)^2
;
;NB: One iteration is scaled by the number of input parameters.

function mcmc_fit, my_funct, params, step_size, num_iter, data=data, errors=errors, cov_mat=cov_mat, $
 chi2_hist=chi2_hist, stat_hist=stat_hist, best_params=best_params,  saved_params = saved_params,  $
	  noplot = noplot,  diag_freq = diag_freq,  save_freq = save_freq

if not keyword_set(diag_freq) then diag_freq = -1
diag_freq =  long(diag_freq)
nprint = 5000 ;!!! Temp only.
num_iter = long(num_iter)
nprint = long(nprint)
init_step_size =  step_size
n_params = long(n_elements(params))
if not keyword_set(save_freq) then save_freq = n_params
print, 'Save frequency: ', save_freq
saved_params = fltarr(n_params,num_iter) 
chi2_hist = fltarr(num_iter)
stat_hist = fltarr(num_iter) 
chi2 = call_function(my_funct, params, data=data, errors=errors, chi2_only=chi2_only)
mobility =  fltarr(n_params)
mobtarg =  0.25
mobility[*] = mobtarg
directions = identity(n_params)
cov_mat = fltarr(n_params,n_params)

for i=long(0),num_iter*long(save_freq)-1 do begin
 chi2_hist[i/save_freq] = chi2_only
 stat_hist[i/save_freq] = chi2
 saved_params[*,i/save_freq] = params
 if (randomu(seed) lt 0.5) then inc = -1.0 else inc = 1.0
 index = randomu(seed)*n_params
 new_params = params + directions[*, index]*step_size[index]*inc
 new_chi2 = call_function(my_funct, new_params, data=data, errors=errors, chi2_only=new_chi2_only)
 if (randomu(seed) lt exp((chi2-new_chi2)/2.0)) then begin
   params = new_params
   chi2 = new_chi2
   chi2_only = new_chi2_only
   mobility[index] =  0.99*mobility[index] + 0.01
 endif else mobility[index] =  0.99*mobility[index]
 if ((i+1) mod 50 eq 0) then begin
  step_size *= (1+ 0.5*(mobility-mobtarg)/mobtarg) 
  if ((diag_freq eq -1) or (i+1) lt diag_freq*save_freq) then step_size =  step_size < 10.*init_step_size
;  print,  step_size ;Lets get this empirically.
 endif
 if ((i+1) mod (nprint*save_freq) eq 0) then begin
   print,  'Done iteration: ', (i+1)/save_freq
   print,  params,  format = '('+strtrim(save_freq, 2)+'F11.2)'
   print,  mobility
   if not keyword_set(noplot) then plot,  chi2_hist,  yr = [0, max([100, chi2_hist[i/save_freq-500]])]
 endif
 if (diag_freq ne -1) then if ((i+1) mod (diag_freq*save_freq) eq 0) then begin
   print,  'Diagonalizing covariance...'
   print,  step_size,  format = '('+strtrim(n_params, 2)+'F11.4)'
   current_ix =  (i+1)/save_freq
   recent_params = saved_params[*, current_ix-diag_freq:current_ix-1]
   mn_params = total(recent_params,2)/diag_freq
   mn_params = rebin(mn_params, n_params, diag_freq)
   cov_mat = (mn_params-recent_params)#transpose(mn_params-recent_params)
   cov_mat = cov_mat/diag_freq
   step_size = la_eigenql(cov_mat,  eigenvectors = directions)
   step_size =  step_size > 1e-6 ;Temp line!!!
   step_size =  0.5*sqrt(step_size) ;Required for curvey minima.
   print,  step_size,  format = '('+strtrim(n_params, 2)+'F11.4)'
   print,  'Continuing...'
 endif
endfor

;Now find the best parameters...
w = where(chi2_hist eq min(chi2_hist))
best_params = saved_params[*,w[0]]

mn_params = total(saved_params,2)/num_iter
for i = long(0),num_iter-1 do cov_mat = cov_mat + (mn_params-saved_params[*,i])#transpose(mn_params-saved_params[*,i])
cov_mat = cov_mat/num_iter

return, mn_params
end
