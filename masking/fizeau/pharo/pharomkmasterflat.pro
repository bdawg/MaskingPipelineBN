;; ----------------------------------------------------------------------------
;;                  make a pharo master flat field
;; ----------------------------------------------------------------------------

pro pharomkmasterflat, filenums, ddir=ddir, outdir=outdir, prefix=prefix, $
                       suffix=suffix, refbias=refbias

if NOT keyword_set(ddir) then ddir = './'
if NOT keyword_set(outdir) then outdir = './'
if NOT keyword_set(prefix) then prefix = 'ph'
if NOT keyword_set(suffix) then suffix = '.fits.gz'

nfiles = n_elements(filenums)
;; check the mode and dimensions of the first file
n = filenums[0]      
fname=prefix+string(n,format="(I4.4)")+suffix

a = readfits(ddir + fname , head0)
mode = pharomode(head0)
print, 'WRITEMOD = ' + mode
naxis1 = 2 * sxpar(head0, 'NAXIS1')
naxis2 = 2 * sxpar(head0, 'NAXIS2')

flats = fltarr(naxis1, naxis2, nfiles)

for i=0,nfiles-1 do begin
  ;; ---------------------- read each file ------------------------------
  fname=prefix+string(filenums[i],format="(I4.4)")+suffix
  print, "reading : "+fname
  a = pharoreadfits(ddir + fname, head, /specklecube, refbias = refbias)

  ;; ---------- if NPAUSEFR != 0, fowler sampling was used -------------
  if ((size(a))[0] gt 2) then begin
    zz = sxpar(head, 'NAXIS3')
    zz = zz/4 -1
  endif else begin
    zz = 1
  endelse

  ;; ------- check that those files can actually be combined ------------
  thismode=pharomode(head)
  print, 'WRITEMOD = ' + thismode
  if thismode ne mode then message,'PHAROMKMASTERFLAT: modes do not match'

  ;; ------- for each pairwise substraction, divide by its median -------
  for j = 0, zz-1 do begin
    medi = median(a[*, *, j])
    a[*, *, j] = a[*, *, j] / medi
  endfor
  
  ;; --------------- average all pairwise substractions -----------------
  ;;    this produces a flat for each filed opened
  for ix = 0, naxis1-1 do begin
    for iy = 0, naxis2-1 do begin
      flats[ix,iy,i] = median(a[ix,iy,*]) ;; median instead of mean ??
    endfor
  endfor
endfor

;; ------------------------ average all flats ---------------------------
masterflat = fltarr(naxis1, naxis2) ;; the flat
sigmaflat  = fltarr(naxis1, naxis2) ;; stability of pixels
for ix = 0, naxis1-1 do begin
  for iy = 0, naxis2-1 do begin
    masterflat[ix, iy] = median(flats[ix, iy,*])
    sigmaflat[ix, iy]  = stddev(flats[ix, iy,*])
  endfor
endfor

;; ----- make a mask to filter some of the bad pixels ------
mask = 1.0 - fltarr(naxis1, naxis2)
;; set here pixels to keep. ex: mask[0,0] = 0
;;mask[465:510,395:440] = 0.0 ;; the black spot
;;mask[*,980:*] = 0.0       ;; the black bar

;; ------------------ bad pixels map --------------------
info = file_info(outdir+'badpix.fits')

if info.exists then begin ;; ---- bad pixel map already exists: open it
  print, 'BAD PIX MAP ALREADY EXISTS'
  pixmap = readfits(outdir+'badpix.fits')

endif else begin ;; ------------- bad pixel map doesn't exist: create it
  pixmap = fltarr(naxis1, naxis2)
  pixmap[*] = 1.0
  print, 'BAD PIX MAP WILL BE CREATED'
endelse

num0 = total(1-pixmap)
;; ------ 1. flag pixels with |gain| > 4.0 * average |gain| ------
medi      = median(masterflat[where(mask gt 0.0)])      ;; exclude hot corners
threshold = 4.0 * stdev(masterflat[where(mask gt 0.0)]) ;; from the stats
temp = where(abs(masterflat-medi) gt threshold)
if temp[0] ne -1 then begin
  pixmap[where(abs(masterflat-medi) gt threshold)] -= 1.0
endif
print, '>>>>> BAD GAIN: ', total(1-pixmap)-num0, ' new pixels'

;; ------ 2. flag pixels with fluctuation > 5.0 * average fluctuation ------
threshold = 5.0 * median(sigmaflat[where(mask gt 0.0)]) ;; exclude weird stuff
temp = where(sigmaflat gt threshold)
if temp[0] ne -1 then begin
  pixmap[where(sigmaflat gt threshold)] -= 1.0
endif
print, '>>>>> UNSTABLE GAIN: ', total(1-pixmap)-num0, ' new pixels'

;; ------ save the map -----
temp = where(mask lt 1.0)
if temp[0] ne -1 then begin
  pixmap[where(mask lt 1.0)] = 1.0
endif
writefits, outdir + 'badpix.fits', pixmap, head0
print, 'BAD PIX MAP UPDATED: ', total(1-pixmap)-num0, ' new pixels added'

;; -------- add keywords to the fits file header before saving ----------
sxaddhist, 'PHAROMKMASTERFLAT.PRO: '+systime()+userid(), head0
sxaddhist, '$Id: pharomkmasterflat.pro,v 1.8 2007/06/13 00:13:18 mireland Exp $',$
           head0
sxaddhist, 'ddir='+ddir, head0
sxaddhist, 'filenums='+strc(filenums), head0
color = sxpar(head,'FILTER')
writefits, outdir+'flatField.'+strc(color)+'.fits', masterflat, head0
end
