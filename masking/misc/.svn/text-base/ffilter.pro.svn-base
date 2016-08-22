;function that filters data in the Fourier plane

function ffilter, data, remove=remove

if (keyword_set(remove) eq 0) then remove = 0

data_fft = fft(data,-1)

plot, modsq(data_fft), xrange=[5,n_elements(data)/2], xstyle=1;, yrange=[0,0.0002], xstyle=1
print, 'Click on left then right of range to filter.'

cursor, x1, y
wait, 0.2
cursor, x2, y
wait, 0.2
x1 = floor(x1)
x2 = floor(x2)

if remove then begin
 data_fft[x1:x2] = 0.
 data_fft[n_elements(data)-x2:n_elements(data)-x1] = 0.
endif else begin
 data_fft[0:x1] = 0.
 data_fft[n_elements(data)-x1:n_elements(data)-1] = 0.
 data_fft[x2:n_elements(data)-x2] = 0.
endelse

return, float(fft(data_fft,1))
end
