;; ACC's 2012 edit of PGT etc's mdark.pro . This version suggests
;; upper and lower cutoffs based on the number of standard deviations
;; from the median counts. This could in theory be completely
;; automated, but it is probably best to check the limits each
;; time.
;;
;; Note that the default cutoff set by ACC is 5 standard deviations.
;; This was chosen for no particular reason other than it looked about
;; what I was normally using and can be changed with the sigma keyword.
;; And it obviously assumes it is roughly Gaussian
;;
;; Setting full auto will not even prompt you to check the limits. Use
;; at your own risk!
;;
;; This program has become a general use program for making flats and 
;; finding bad pixels

pro mdark2,files,flat,bad,cutoff,sigma=sigma,full_auto=full_auto

if (keyword_set(files) eq 0) then begin
 print,'USAGE:'
 print,'mdark2,files,flat,bad,cutoff'
 return
endif

nflats=n_elements(files)
if not keyword_set(sigma) then sigma=5

for i=0,nflats-1 do begin
   filename=files(i)
   nflat=float(reform(readfits(filename,head,/silent)))
   sz=size(nflat)
   if (i eq 0) then begin
       flat=fltarr(sz[1],sz[2]) 
       dumbad=flat
   endif

   ;;Anthony's quick way
   hist=histogram(nflat,locations=xaxis)
   plothist,nflat,xhist,yhist,/noplot
   med=median(nflat)
   sd=sqrt(variance(nflat))
   
   ;;truncate the histogram at a certain number of standard deviations
   ;;from the median
   max_ok=med+sigma*sd
   min_ok=med-sigma*sd
   
  
   ;;Now display a plot so we can check that it worked
   plot,xaxis,hist,title='Click outside area to accept bounds, or inside to switch to manual bounds'
   oplot,[min_ok,min_ok],[0,5*max(hist)],line=1
   oplot,[max_ok,max_ok],[0,5*max(hist)],line=2
   if keyword_set(full_auto) then v1 =0 else cursor,v1,a,/device
   wait,0.2
   
   ;;the edge of the plot should be 60 in device coordinates, so use
   ;;that
   if v1 lt 60 then good=where(nflat lt max_ok and nflat gt min_ok ,complement=w) else begin
       try_again:
       plothist,nflat,tit='Click on left then right Edge of allowed range'
       cursor,v1,a &wait,.6
       print,' Lower Bound ',v1
       cursor,v2,a & wait,.6
       print,' Upper Bound ',v2
       good=where(nflat gt v1 and nflat lt v2)
       if good[0] eq -1 then begin
           print,'You selected zero pixels! Try again!'
           goto,try_again
       endif
       plothist,nflat[good],tit='SECOND CHANCE: Click on left then right Edge of allowed range'
       cursor,v1,a &wait,.6
       print,' Lower Bound 2',v1
       cursor,v2,a & wait,.6
       print,' Upper Bound 2',v2
       good=where(nflat gt v1 and nflat lt v2)
       if good[0] eq -1 then begin
           print,'You selected zero pixels! Try again!'
           goto,try_again
       endif
       w=where(nflat lt v1 or nflat gt v2)
       plothist,nflat[good],tit='Final Histogram'
   endelse
   flat=flat+nflat
   dumbad(w)=dumbad(w)+1
endfor

bad=where(dumbad gt cutoff)

flat=flat/nflats

print,n_elements(bad),' bad pix flagged'

end
