;+
; sky_subtract_nirc,image,radius,cent=cent,trim=trim,stats=stats,plotit=plotit,
;                   printit=printit,bin=bin,bad_mask=bad_mask
; 
;  This will take an image, mask out a circular portion of it, then fit a 
;  gaussian to whats left.
;  Then it will subtract off the position of the maximum of the gaussian.
;  It will optionally Output what was subtracted.
;  bad_mask is a mask of bad pixels (0 means good, 1 means bad) 
;
;  STATS: will contain (area of gaussian, peak position, sigma)
;
;  PGT - fixed bad bug which made convergence screw up for bright stars.
;        problem was that initial estimates given to curvefit came from
;        whole frame rather than just edges .....
;  
;  PGT - added subtr keyword, to give the option not to subtract the sky level
;  sky_subtract_nirc for new software release. Removed leak keywords  PGT 05Dec03
;  PGT - trim keyword cuts a fraction trim off the edge (which are bright for nirc2) 
;  PGT 30Apr04 - If binsize not set, use adaptive bin

pro sky_subtract_nirc,image1,radius,cent=cent,trim=trim,plotit=plotit,printit=printit,$
	binsize=binsize, stats=stats,bad_mask=bad_mask,subtr=subtr
;	

if (keyword_set(trim) ne 0) then begin
   dimx=(size(image1))[1] & dimy=(size(image1))[2]
   newdimxy=min([dimx,dimy]*trim)
   testim=grabnxn(image1,newdimxy,x=dimx/2,y=dimy/2)
   radius=radius*trim
endif else  testim=image1

info=size(testim)
x=info(1)
y=info(2)
if (keyword_set(bad_mask) eq 0) then bad_mask=fltarr(x,y)
if (keyword_set(subtr) eq 0) then subtr=0 ; subtract sky by default

if (keyword_set(cent) eq 0) then begin
  locate_peak,smooth(testim*(1-bad_mask),10),xpeak,ypeak
endif else begin
  xpeak=cent(0)
  ypeak=cent(1)
endelse

r=radius
image=testim
; imagenoise=sky_noise(image,imagelevel)  **** PGT Bugfix removed this line
;                                       will not converge when star is bright
mask=fltarr(x,y)
make_circle,mask,x/2,y/2,radius
mask=shift(mask,xpeak-x/2,ypeak-y/2)
index=where(mask eq 0 and bad_mask eq 0)
if (keyword_set(binsize) eq 0) then binsize=stdev(image(index))/20.
plothist,image(index),x1,v2,bin=binsize,/noplot
imagenoise=sky_noise(image(index),imagelevel)  ;**** PGT Bugfix. Moved line to here ...

b2=[total(v2)*binsize,imagelevel,imagenoise]
weights=replicate(1.0,n_elements(v2))
parfit2=curvefit(x1,v2,weights,b2,sig2,function_name='skygauss')

if (keyword_set(printit) eq 1) then begin 
  print,'For Ima: Area,Disp,Sigma:',b2,'(',sig2,')'
endif
if (keyword_set(plotit) eq 1) then begin
  plot,x1,v2,tit="Sky Subtraction",xr=[b2(1)-3*b2(2),b2(1)+3*b2(2)]
  oplot,x1,parfit2,line=2
endif

stats=b2
if (subtr lt 1.0) then begin    
  image1=image1-stats(1)
  if (keyword_set(plotit) eq 1) then print,'Performing Sky Subtraction'
endif

end
	
