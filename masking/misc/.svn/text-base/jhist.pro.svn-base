function jhist,im1a,nz=nz
if (keyword_set(im1a) eq 0) then begin
  print,' image= jhist(image,nz=nz) '
 return,-1 
endif 
im1=im1a
if (keyword_set(nz) eq 1) then goto,skip

; original by JDM in 1998
; 25Feb99 JDM  Added nz=33 feature to allow
; the dynamic rnage not to wasted by NOISE.



n=n_elements(im1)
random_order=sort(randomn(seed,n))
im=im1
im(*)=im1(random_order)

in=sort(im)
n=n_elements(im)


col=findgen(n) & col=col*255/max(col)
col=fix(col)

newim=im
newim(random_order(in))=col

return,newim

skip:  ; Go here for nz treatment

result=histogram(im1a,binsize=nz,reverse_indices=ri)
in=where(result ne 0,ct)
cols=findgen(ct)*255/(ct-1)
for i=0l,ct-1l do begin
   a=ri(in(i))
   b=(ri( (in(i))+1)-1)
;   if (a ne b) then $
      im1(ri(a:b)) = cols(i) 
;	else $
;      im1(ri(a)) = cols(i)
;   print,i
     
endfor
return,im1

end


