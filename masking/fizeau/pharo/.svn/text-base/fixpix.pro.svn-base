pro fixpix,imgs,badpix,outimgs, $
	   npix=npix,weight=weight, $
	   noise=noise,sigma=sigma,dc=dc, $
	   silent=silent

;+
; given a image or stack of images and a bad pixel
; mask, will fill in bad pixels by finding the NPIX nearest
; good pixels, toss the highest and lowest ones of the bunch,
; and then arithmatically average. 
;
; bad pixel list is processed in order array is stored in
; memory, i.e. row by row, and is defined by:
;	 0 = bad pixel
;    not 0 = good pixel
;
; NPIX = # of adjacent pixels used for correction (default = 8)
; /weight: weight adjacent pixels by inverse of their distances
;	   in averaging process
;
; checks to see if some pixels are equidistant; if so,
; use all of them, even if the total # of pixels then exceeds NPIX
;
; WARNINGS:
;  - will work for entire bad columns/rows, but 
;    may not be the most sensible thing to use
;  - do not pass it the same array for input & output as this
;    might corrupt the correction algorithm
;
; 7/3/95 MCL
;
; added /noise: replace bad pixels with random noise with
; user defined sigma and dc level
; 9/24/95 MCL
;
; badpix can now be a 2d or 3d array
; 4/18/96 MCL
;
; uses MESSAGE now instead of regular PRINT
; 04/25/01 MCL
;-


if keyword_set(npix) eq 0 then npix = 8
if n_params() lt 2 then begin
	print,'pro fixpix,imgs,badpix,outimgs,'
	print,'           [npix=',strc(npix),'],[weight],'
	print,'           [noise],[sigma=],[dc=],[silent]'
	retall
endif

;if n_elements(outimgs) ne 0 then $
;  if (total(outimgs-imgs) eq 0) then $
;	message,'original and output images are the same!'
outimgs = imgs		; don't change originals!

if keyword_set(noise) and not(keyword_set(sigma)) then begin
	read,'enter sigma for random noise: ',sigma  
	read,'enter dc level for noise: ',dc
endif
if keyword_set(noise) and not(keyword_set(dc)) then dc=0

sz = size(imgs)
if sz(0) eq 2 then begin
	nimg = 1 
	img = imgs
endif else if sz(0) eq 3 then begin
	nimg = sz(3) 
	img = imgs(*,*,0)
endif else begin
	print,'** input image is not 2 or 3-d! **'
	return
endelse
szb = size(badpix)
if szb(0) eq 3 and szb(3) ne sz(3) then begin
    message, '3d badpix file different size than image file', $
      /info 
    return
endif


; loop through images
for j=0,nimg-1 do begin

    if szb(0) eq 3 then bp = badpix(*, *, j) $
    else bp = badpix
    wbp = where(bp eq 0.0, nbad)
    if keyword_set(silent) eq 0 then begin
        if (nimg gt 1) then print, format = '("image ",A,": ")', strc(j)
;        print, format = '($,"image:" )'
;	print,format='(" ",A)',strc(j)
        message, strc(nbad)+' bad pixels to fix', /info
  endif
  if sz(0) eq 3 then img = imgs(*,*,j)

  imoffset = j * sz(1) * sz(2)

; loop through pixels & replace with noise
  if keyword_set(noise) then begin

    for i=0L,nbad-1 do begin
	newval = randomn(seed) * sigma + dc
	outimgs(wbp(i)+imoffset) = newval
    endfor

; or interpolate
  endif else begin
	
    for i=0L,nbad-1 do begin

;	default search box is 2*dd+1 on a side
	dd = 2
	found = -1

;	determine search region
	repeat begin
		y = floor(wbp(i)/sz(1))
		x = wbp(i) - y*sz(1)
		x0 = x-dd > 0
		x1 = x+dd < (sz(1)-1)
		y0 = y-dd > 0
		y1 = y+dd < (sz(2)-1)
		bpcut = bp(x0:x1,y0:y1)
		cut = img(x0:x1,y0:y1)
		wgood = where(bpcut ne 0.0,ngood)
		if ngood lt npix then begin
;			print,'expanding dd'
			dd = dd + 2
		endif else found=1
	endrep until (found ne -1)

;	calculate distances to adjacent good pixels 
	dist_circle,distarr,2*dd+1,dd,dd
	gdist = distarr(wgood)
	gpix = cut(wgood)

;	sort good pixels by distance
	ss = sort(gdist)
	gdist = gdist(ss)
	gpix = gpix(ss)


;	accounting for pixels with the same distance at the edge
	mm = where(gdist(npix-1:*) eq gdist(npix-1),nn)
	nn = nn - 1
;	if nn gt 0 then print,'  multiplicty = ',strc(nn+1)

;	error checking
;	if i eq 0 then begin 
;		for k=0,npix+nn-1 do print,gdist(k),gpix(k)
;		print,x,y
;		print,cut
;	endif

;	get the values of the relevant pixels
	gpix = gpix(0:npix+nn-1)
	ss2 = sort(gpix)
	gpix = gpix(sort(gpix))
	gdist = gdist(ss2)

;	calculate new pixel value, tossing the highest
;	and lowest pixels of the bunch, then weighting
;	by the inverse of the distances if desired
	if keyword_set(weight) then begin
	  newval = total(gpix(1:npix+nn-2)/gdist(1:npix+nn-2))
	  newval = newval / total(1./gdist(1:npix+nn-2))
	endif else begin
	  newval = total(gpix(1:npix+nn-2)) / (npix+nn-2.0)
	endelse

;	more error checking
	if i eq 0 and keyword_set(silent) eq 0 then $
          print, '  oldval=', strc(outimgs(wbp(i)+imoffset)), ' ' + $
          'newval = ', strc(newval) 

	outimgs(wbp(i)+imoffset) = newval

    endfor
 
  endelse

endfor
	


end
