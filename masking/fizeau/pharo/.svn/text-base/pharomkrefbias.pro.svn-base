; make a pharo dark from a sequence of frames
;$Id: pharomkrefbias.pro,v 1.6 2007/07/19 21:47:21 mireland Exp $
;$Log: pharomkrefbias.pro,v $
;Revision 1.6  2007/07/19 21:47:21  mireland
;Unsure why mkrefbias has changed. Fixed a HUGE bug in the mrg file
;creation from calibrate_v2_cp.pro
;
;Revision 1.5  2007/06/13 00:13:18  mireland
;Added a capability to only call pixels bad if they are bad in multiple frames.
;
;Revision 1.3  2006/05/04 18:45:01  frantzm
;now updates the badpix map with dead/hot and unstable pixels
;
;Revision 1.2  2006/02/09 22:23:25  frantzm
;pharoreadfits now uses a correct bad pixel treatement procedure
;
;Revision 1.1.1.1  2005/12/19 05:15:05  mireland
;Initial release of masking code.
;
;Revision 1.3  2003/12/19 23:03:40  jpl
;*** empty log message ***
;
;Revision 1.2  2003/11/07 19:55:47  jpl
;*** empty log message ***
;

pro pharomkrefbias, filenums, ddir=ddir, outdir=outdir, prefix=prefix, $
                    suffix=suffix
   

if NOT keyword_set(ddir) then ddir = './'
if NOT keyword_set(outdir) then outdir = './'
if NOT keyword_set(prefix) then prefix = 'ph'
if NOT keyword_set(suffix) then suffix = '.fits.gz'

nfiles = n_elements(filenums)

; check the mode of the first file
n = filenums[0]      
fname=prefix+string(n,format="(I4.4)")+suffix
print, 'image: ' + ddir + fname
a = readfits(ddir + fname ,head0)
mode = pharomode(head0)
print, mode
naxis1 = sxpar(head0,'NAXIS1')
naxis2 = sxpar(head0,'NAXIS2')
naxis3 = sxpar(head0,'NAXIS3')

refbiases = fltarr(naxis1,naxis2,naxis3,nfiles) 

for i=1,nfiles-1 do begin  
  fname=prefix+string(filenums[i],format="(I4.4)")+suffix
  print, "reading :"+fname
  a = readfits(ddir + fname, head)
  
  thismode=pharomode(head)
  print, thismode
  if thismode ne mode then message,'PHAROMKREFBIAS: modes do not match'
  refbiases[*,*,*,i] = a
endfor


if (nfiles gt 1) then begin
  refbias = fltarr(naxis1,naxis2,naxis3)      ;; the master dark
  sigmadark = fltarr(naxis1, naxis2, naxis3)  ;; stability of pixels
  ;; ------------------------- median algorithm -------------------------------

  ti = systime(1)
  for i=0,naxis3-1 do begin
    print, string(i+1,format="(I2)"), '/', string(naxis3,format="(I2)"), $
           ': time past ', string(systime(1)-ti, format="(F6.2)"), ' seconds'
    medarr, reform(refbiases[*,*,i,*]), r
    refbias[*,*,i] = r

    for ix = 0, naxis1-1 do begin
      for iy = 0, naxis2-1 do begin
        sigmadark[ix, iy, i] = stddev(refbiases[ix, iy, i, *])
      endfor
    endfor
  endfor
  sigmadark = sigmadark/median(sigmadark)
endif else begin
  refbias = reform(refbiases) ;; make a dark with just 1 file
  sigmadark[*,*,*] = 0.0      ;; sounds suicidal no ?
endelse

;; ---------------------- mask the last column ----------------------------
;; even with 20 darks, the last column of each median quadrant is full
;; of fake bad pixels, which therefore biases the statistics of the frame.
;; Here, we have to cheat a bit ...

for i=0,naxis3-1 do begin
  medi = median(refbias[*,*,i])
  refbias[naxis1-1,*,i] = medi   ;; the last column is excluded from bad piels
  sigmadark[naxis1-1,*,i] = 0.0 ;; idem here ...
endfor

;; -------------------------- bad pixels map -------------------------------
info = file_info(outdir + 'badpix.fits')

if info.exists then begin ;; ---- bad pixel map already exists: open it
  pixmap = readfits(outdir + 'badpix.fits')
  print, 'BAD PIXEL MAP ALREADY EXISTS'
endif else begin ;; ------------- bad pixel map doesn't exist: create it
  pixmap = fltarr(1024, 1024)
  pixmap[*] = 1.0
  print, 'BAD PIXEL MAP WILL BE CREATED'
endelse

num0 = total(1-pixmap)
;; ---------------- update the bad pixel map and save --------------------
nlayer = naxis3/4
px = 512
py = 512

minimap    = fltarr(2*naxis1, 2*naxis2) ;; prepares a mini badpix map
dummy      = fltarr(1024, 1024)         ;; and a full size temp map

for i = 0 , nlayer-1 do begin
  minimap[*] = 1.0
  layer      = pharojoin(refbias[*,*,4*i:4*i+3])
  sigmaLayer = pharojoin(sigmadark[*,*,4*i:4*i+3])

  ;; ----- 1. flag pixels with |value| > 4.0 * average |value| -----
  medi      = median(refbias[0:naxis1-2,*,4*i:4*i+3]) ;; last column excluded
  threshold = 4*median(abs(refbias[0:naxis1-2,*,4*i:4*i+3])) ;; from the statistics
  w = where(abs(layer - medi) gt threshold)
  if (w[0] ne -1) then minimap[w] = 0.0
 
  ;; ----- 2. flag pixels with fluctuation > 4.0 * average fluctuation -----
  w = where(sigmaLayer gt (4.0 >  1 + 10./sqrt(nfiles)))
  if (w[0] ne -1) then minimap[w] = 0.0

  ;; projects the minimap on the real map
  dummy[*]  = 1.0               ;; init the dummy map

  minis = pharosplit(minimap)
  dummy[px:px+naxis1-1, 0:naxis2-1]     = minis[*,*,0]
  dummy[0:naxis1-1,     0:naxis2-1]     = minis[*,*,1]
  dummy[0:naxis1-1,     py:py+naxis2-1] = minis[*,*,2]
  dummy[px:px+naxis1-1, py:py+naxis2-1] = minis[*,*,3]
  minis = 0                     ;; free memory
  pixmap += dummy-1       ;; add new bad pixels to the map
endfor

writefits, outdir + 'badpix.fits', pixmap
print, 'BAD PIXEL MAP UPDATED: ', string(total(1-pixmap)-num0,format="(I5)"),$
       ' new pixels added'

minimap = 0     ;; free memory
dummy   = 0     ;; free memory
pixmap  = 0     ;; free memory

;; --------- add history to the fits file header before saving -----------
print, 'saving reference bias'
sxaddhist, 'PHAROMKREFBIAS.PRO: '+systime()+userid(), head0
sxaddhist, '$Id: pharomkrefbias.pro,v 1.6 2007/07/19 21:47:21 mireland Exp $',head0
sxaddhist, 'ddir='+ddir, head0
sxaddhist, 'filenums='+strc(filenums), head0
writefits, outdir+pharorefbiasname(head0), refbias, head0
end
