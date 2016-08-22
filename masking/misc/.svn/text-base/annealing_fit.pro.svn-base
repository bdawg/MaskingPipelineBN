;This fuction returns a mean and covariance matrix of a bunch of
;parameters, given some input data. Two variables, data and errors can
;be passed through this routine to the chi^2. It is of course
;important that 'data' and 'errors' are both passed by reference, ie
;without the use of subscripts...
;
;my_funct:     A string containing the name of the chi^2 calculation
;              function.
;params:       Initial parameters
;step_size:    Step size at a temperature of 1. Scaled by sqrt(temp) 
;start_temp:   For simulated annlealing-style calculations, start_temp
;              divides the initial chi^2 by this number
;cooling_time: Number of iterations for 1/e decrease in
;              temperature. Cooling stops altogether after 5 times
;              this time. 
;num_iter:     Number of iterations after cooling has completed.
;NB: One iteration is scaled by the number of input parameters.

function annealing_fit, my_funct, params, step_size, start_temp, cooling_time, num_iter, data=data, errors=errors, cov_mat=cov_mat, fixed=fixed

n_params = n_elements(params)
if keyword_set(fixed) then fixed = intarr(n_params)
saved_params = fltarr(n_params,num_iter) 
chi2_hist = fltarr(num_iter + cooling_time*8)
temp = start_temp
chi2 = call_function(my_funct, params, data=data, errors=errors)/(1.0+temp)

for i=long(0),(num_iter + cooling_time*5)*long(n_params)-1 do begin
 chi2_hist[i/long(n_params)] = chi2
 oldtemp = temp
 if (i lt 5*cooling_time*long(n_params)) then begin
       temp = start_temp*exp(-i/float(n_params)/cooling_time) + 1.0
       ;chi2 = chi2*oldtemp/temp
     endif else begin
     temp = 1.0 ;ie temperature, not necessarily a `temporary' variable
     if (i mod long(n_params) eq 0) then saved_params[*,i/long(n_params) - 5*cooling_time] = params
 endelse
;This next line should not be needed
chi2 = call_function(my_funct, params, data=data, errors=errors)/temp
 if (randomu(seed) lt 0.5) then inc = -1.0 else inc = 1.0
 index = randomu(seed)*n_params
 new_params = params
 new_params[index] = new_params[index] + step_size[index]*inc
 new_chi2 = call_function(my_funct, new_params, data=data, errors=errors)/temp
 if (randomu(seed) lt exp(chi2-new_chi2)) then begin
   params = new_params
   chi2 = new_chi2
 endif
endfor

mn_params = total(saved_params,2)/num_iter
cov_mat = fltarr(n_params,n_params)
for i = 0,num_iter-1 do cov_mat = cov_mat + (mn_params-saved_params[*,i])#transpose(mn_params-saved_params[*,i])
cov_mat = cov_mat/num_iter

;Temp stuff specific to fitting a binary orbit...
orbit_fit = fltarr(2,1000)
for i = 0,999 do orbit_fit[*,i] = binary_position([mn_params[0]+i*mn_params[1]/1000.0,mn_params[1:6]])
ofit = fltarr(2,24)
for i = 0,(size(data))[2]-1 do ofit[*,i] = binary_position([data[0,i]-mn_params[0],mn_params[1:6]])
plot, orbit_fit[0,*]*cos(orbit_fit[1,*]*!pi/180), orbit_fit[0,*]*sin(orbit_fit[1,*]*!pi/180)
oplot, data[1,*]*cos(data[2,*]*!pi/180), data[1,*]*sin(data[2,*]*!pi/180), psym=5
oplot, ofit[0,*]*cos(ofit[1,*]*!pi/180), ofit[0,*]*sin(ofit[1,*]*!pi/180), psym=4

stop

dummy = 0
return, mn_params
end
