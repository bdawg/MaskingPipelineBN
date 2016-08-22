function rebinmed, img, xnew,ynew, reducf=reducf
    sz=size(img)
    xold=sz[1]
    yold=sz[2]

    if not(keyword_set(reducf)) then reducf=xold/xnew
 ;   print, "reducf = "+strc(reducf)
    out=fltarr(xnew,ynew)
    for j=0,ynew-1 do begin
        for i=0,xnew-1 do begin
            x0 = i*reducf
            x1 = (i+1)*reducf -1
            y0 = j*reducf
            y1 = (j+1)*reducf -1

;            print, "newpix/"+strc(x0)+":"+strc(x1)+","+strc(y0)+":"+strc(y1)+" -> "+strc(i)+","+strc(j)
            out[i,j] = median( img[ x0:x1, y0:y1 ] )
;            print, total(img[ x0:x1, y0:y1 ] ), out[i,j]
       endfor
    endfor
    return,out
end

function phfindstarquad, thedata, reducf=reducf

    sz=size(thedata)

    nframes = sz[3]

    if not(keyword_set(reducf)) then reducf=min(sz[1:2])/32
        
    sz = sz/reducf
    LA = sz * reducf -1   ;may have to drop edges of images.

    imtmp = fltarr(sz[1],sz[2],nframes,4)
    for i=0, nframes-1 do begin
        for j=0,3 do $
            imtmp[*,*,i,j] = rebinmed( thedata[ 0:LA[1], 0:LA[2], i,j ], sz[1], sz[2]  )
    endfor
    fluxes = total(total(total(imtmp,1),1),1)
    maxflux = max(fluxes,i)
    return, i
end
