;;Author ACC 2013
;; Makes flats from twilight files. This program asks you
;; what files to use as on and off, then calls mdark2.pro
;; to do the actual flat making
pro twiflats2,twi_files,flat,bad,cutoff,sigma=sigma,full_auto=full_auto

;;first, we need to ask which files are "on" and which are "off" files
im=readfits(twi_files[0],/silent)
if (size(im))[0] gt 2 then im=mean(im,dim=3)
twis=im
flux=median(im)
for i=1,n_elements(twi_files)-1 do begin
   im=readfits(twi_files[i],/silent)
   if (size(im))[0] gt 2 then im=mean(im,dim=3)
   twis=[[[twis]],[[im]]]
   flux=[flux,median(im)]
endfor

;;plot median flux and ask which ones to use as on files.
;; In theory you could also automate this (or make a suggestion)
try_again:
plot,flux,title='Select bounds for "on" files',xtitle='File number',ytitle='Median flux'
cursor,x,y,/data & wait,0.2
if x lt 0 then x=0.1
oplot,[x,x],[min([flux,0]),max(flux)],linestyle=1

cursor,x2,y2,/data & wait,0.2
oplot,[x2,x2],[min([flux,0]),max(flux)],linestyle=2
if x2 gt n_elements(flux)-1 then x2=n_elements(flux)-1

if x2 lt x then begin
   print,'You selected a negative range!'
   goto, try_again
endif

;;plot median flux and ask which ones to use as off files.
;; In theory you could also automate this (or make a suggestion)
try_again3:
plot,flux,title='Select bounds for "off" files',xtitle='File number',ytitle='Median flux'
oplot,[x,x],[min([flux,0]),max(flux)],linestyle=1,color=100
oplot,[x2,x2],[min([flux,0]),max(flux)],linestyle=2,color=100
wait,0.2
cursor,x3,y3,/data & wait,0.2
if x3 lt 0 then x3=0.1
oplot,[x3,x3],[min([flux,0]),max(flux)],linestyle=1

cursor,x4,y4,/data & wait,0.2
oplot,[x4,x4],[min([flux,0]),max(flux)],linestyle=2
if x4 gt n_elements(flux)-1 then x4=n_elements(flux)-1

if x4 lt x3 then begin
   print,'You selected a negative range!'
   goto, try_again3
endif

on_files=twi_files[floor(x):ceil(x2)]
off_files=twi_files[floor(x3):ceil(x4)]

;;now call mdark2
mdark2,on_files,flat,flbad,cutoff,sigma=sigma,full_auto=full_auto
mdark2,off_files,dark,dkbad,cutoff,sigma=sigma,full_auto=full_auto

big_bad=0*flat
big_bad[flbad]=1
big_bad[dkbad]=1
bad=where(big_bad gt 0,complement=good)
flat[bad]=median(flat[good])
dark[bad]=median(dark[good])
flat=flat-dark
end
