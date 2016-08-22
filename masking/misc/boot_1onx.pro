;This procedure uses a balanced bootstrap to estimate 
;1/<x>, with x an input vector.

function boot_1onx, x, var=var, niter=niter

if (keyword_set(niter) eq 0) then niter = 400
nx = n_elements(x)
index = intarr(nx,niter)
for i = 0,niter-1 do index[*,i]=indgen(nx)
in1 = randomu(seed,nx*niter*4)*nx*niter
in2 = randomu(seed,nx*niter*4)*nx*niter

for j = long(0),niter*nx*4-1 do begin
 temp = index[in1[j]]
 index[in1[j]] = index[in2[j]]
 index[in2[j]] = temp
endfor
results = fltarr(niter)

for i = 0,niter-1 do results[i] = 1/mean(x[index[*,i]])
result_mean = mean(results)
var = total((results-result_mean)^2)/(niter-1.0)
stop

return, result_mean
end
