;This fuction returns a mean and covariance matrix of a bunch of
;parameters, given some input data. The variable data (which can be a
;structure)  is passed to my_funct, and the variable deriv must be returned.
;
; We use the Hamiltonian Markov Chain algorithm here.
;
;my_funct:     A string containing the name of the chi^2 calculation
;              function.
;params:       Initial parameters - must be reasonable within errors.
;step_size:    Step size for each parameter
;num_iter:     Number of iterations. I think to get a 10% error on
;              covariances you need 400*(sig/step_size)^2
;
;NB: One iteration is scaled by the number of input parameters.

function hmcmc_fit, my_funct, params, step_size, num_iter, data=data, cov_mat=cov_mat, chi2_hist=chi2_hist,$
                   best_params=best_params, mobility=mobility,  saved_params = saved_params

init_step_size =  step_size
n_params = n_elements(params)
saved_params = fltarr(n_params,num_iter) 
chi2_hist = fltarr(num_iter)
chi2 = call_function(my_funct, params, data=data, deriv = deriv0)
tstep =  0.5 ;!!! Check this: smaller than (d^2H/dx^2)^(-0.5)

for i=0,num_iter-1 do begin
 ;Sort out  histories...
 chi2_hist[i] = chi2
 saved_params[*,i] = params
;Find the derivative at the start (!!! maybe this isn't needed if I'm clever)
 chi2 = call_function(my_funct, params, data=data, deriv = deriv0)
 ;Make a randomised velocity vector: (NB: A velocity of 1, or 1 per parameter?)!!!
 u = randomn(seed, nparams)*step_size
 nsteps =  1; + randomu(seed)*10 ;Put this bit in once tested !!!
 for j = 0, nsteps-1 do begin
  ;Make the 'Leapfrog' step:
  ueps2 =  u - tstep/2.*deriv0/2.
  new_params =  params + tstep*ueps2
  new_chi2 = call_function(my_funct, params, data=data, deriv = deriv1)
  ueps =  ueps2 - tstep/2.*deriv1/2.
 endfor
 ;Now for Metropolis decision:
 if (randomu(seed) lt exp((chi2-new_chi2)/2. + total((deriv1/stepsize)^2)/2. - total((deriv0/stepsize)^2)/2) then begin
   params = new_params
   chi2 = new_chi2
   deriv0 = deriv1
 endif

;!!!Printouts: for when things are working !!!
; if ((i+1) mod (1000*n_params) eq 0) then begin
;   print,  'Done iteration: ', (i+1)/n_params
;   print,  params
;   print,  mobility
;   plot,  chi2_hist,  yr = [0, max([100, chi2_hist[i/long(n_params)-500]])]
; endif
endfor

;Now find the best parameters...
w = where(chi2_hist eq min(chi2_hist))
best_params = saved_params[*,w[0]]

mn_params = total(saved_params,2)/num_iter
cov_mat = fltarr(n_params,n_params)
for i = long(0),num_iter-1 do cov_mat = cov_mat + (mn_params-saved_params[*,i])#transpose(mn_params-saved_params[*,i])
cov_mat = cov_mat/num_iter

return, mn_params
end
