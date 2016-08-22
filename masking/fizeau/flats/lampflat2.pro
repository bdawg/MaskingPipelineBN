; 20Feb04
; PGT
;
;;Edited ACC 2012 to make it simpler to use, and to handle flats at
;;longer wavelengths (i.e. using dome flats instead of lamp)
;;
;; Note that this is very similar to mdark2.pro, and there is no need
;; for both. (This could potentially call mdark2, like twiflat2 does)
; USAGE:
;;
;;INPUTS
;;    on: filenames for lamp on
;;   off: filenames for lamp off
;;cutoff: number of times a pixel has to be flagged to be counted as bad
;;
;;OUTPUTS
;;  flat: 
;;   bad: the bad pixel map
;;
;;OPTIONS
;;sigma : this is the number of standard deviations away from the mean
;;        that a pixel has to be before it is called bad. Obviously
;;        assumes that it is roughly Gaussian. Default=5
;;nolamp: For longer wavelengths, no lamp is used and the flats are
;;        actually dome flats, which are dark subtracted
;;        instead. (Darks are .idlvar files, and so need to be handled
;;        differently)
;;nodarks: Maybe we don't even need darks for the longer
;;         wavelength. If dark current is negligible, we can just use
;;         the dome flats
;; Setting full auto will not even prompt you to check the limits. Use
;; at your own risk!

pro lampflat2,on,off,flat,bad,cutoff,sigma=sigma,nolamp=nolamp,no_darks=no_darks,full_auto=full_auto,average_off=average_off

if (keyword_set(on) eq 0 and keyword_set(off) eq 0) then begin
    print,'USAGE:'
    print,'lampflat2,on,off,flat,bad,cutoff,sigma=sigma,nolamp=nolamp,no_darks=no_darks,full_auto=full_auto,average_off=average_off'
    return
endif

 if not keyword_set(sigma) then sigma=5

nflats=n_elements(on)

if keyword_set(average_off) then begin
   av_flatoff=float(reform(readfits(off[0],head,/silent)))
   for i=1,n_elements(off)-1 do begin
      av_flatoff+=float(reform(readfits(off[i],head,/silent)))
   endfor
   av_flatoff/=n_elements(off)
endif

for i=0,nflats-1 do begin
    flaton=float(reform(readfits(on[i],head,/silent)))

    if keyword_set(no_darks) then begin
        flatoff=0.*flaton
    endif else if keyword_set(nolamp) then begin
        restore,off
        flatoff=dark
     endif else if keyword_set(average_off) then begin
        flatoff=av_flatoff
     endif else flatoff=float(reform(readfits(off[i],head,/silent)))
    sz=size(flaton)

    nflat=flaton-flatoff
    if (i eq 0) then begin
        flat=fltarr(sz[1],sz[2]) 
        dumbad=flat
    endif
    
    ;;ACC's automatic way, stolen from mdark2:
    plothist,nflat,xhist,yhist,/noplot
    med=median(nflat)
    sd=sqrt(variance(nflat))
    
    ;;truncate the histogram at a certain number of standard deviations
    ;;from the median
    max_ok=med+sigma*sd
    min_ok=med-sigma*sd
    
    ;;Now display a plot so we can check that it worked
    plot,xhist,yhist,title='Click outside area to accept bounds, or inside to switch to manual bounds'
    oplot,[min_ok,min_ok],[0,5*max(yhist)],line=1
    oplot,[max_ok,max_ok],[0,5*max(yhist)],line=1
    if keyword_set(full_auto) then v1=0 else cursor,v1,a,/device
    wait,0.2
    
    ;;note device coordinate 60 is the corner of the plot
    if v1 lt 60 then good=where(nflat lt max_ok and nflat gt min_ok ,complement=w) else begin
        try_again:
        plothist,nflat, tit=string(i)+': Click the left then right edges of the allowed range.'
        cursor, v1,a &wait,.2
        cursor, v2,a &wait,.2
        print, 'Bounds from '+string(v1)+' to '+string(v2)
        good=where(nflat gt v1 and nflat lt v2)
        if good[0] eq -1 then begin
            print,'Whoops! No pixels selected! Try Again!'
            goto,try_again
        endif
        plothist,nflat[good], tit=string(i)+' Second Chance: Click the left then right edges of the allowed range.'
        cursor, v1,a &wait,.2
        cursor, v2,a &wait,.2
        print, 'Bounds from '+string(v1)+' to '+string(v2)
        good=where(nflat gt v1 and nflat lt v2)
        if good[0] eq -1 then begin
            print,'Whoops! No pixels selected! Try Again!'
            goto,try_again
        endif
        w = where(nflat lt v1 or nflat gt v2)
        plothist, nflat[good], tit="Final Histogram" &wait, .2
    endelse
    
    flat=flat+nflat ;;bad pix get removed later
    dumbad(w)=dumbad(w)+1
endfor

bad=where(dumbad gt cutoff)
flat=flat/nflats

print,n_elements(bad),' bad pix flagged'

end






