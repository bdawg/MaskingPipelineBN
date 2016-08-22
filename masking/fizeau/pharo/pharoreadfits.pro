; $Id: pharoreadfits.pro,v 1.11 2009/09/28 23:54:17 snert Exp $
;
; Description:
;
; read a pharo file and reformat it to be a single frame or
; specklecube
; check writemode
; passes all readfits paramaters 
; 
; doesn't yet check NPAUSEFR for a too long exposure in the middle 
; of the stack.  fowler sampling mode does work.
; todo: add history to for refbias
; 
; use a directory of refbiases if set
;
; Quick explanation of pharo specialty keywords:
; All binary switches:
;  nocds          - don't subtract cumulative frames
;  specklecube    - pairwise subtractions of multiple reads (normally set)
;  quadrants      - if set, gives 4 separate arrays. (default is to tile)
;  ignorenpausefr - Set this to override default STOP if pauseframes are present
;  refbias        - dark file (or directory containing dark to search)
;  flatf          - Flat field (fits)
;  badpixf        - Bad Pix file (fits) 0=bad, non0=good
;  fixpix_rs      - run automated bad pixel fixer.
;  firstbad       - first 'firstbad' frames in each cube are bad
;  display        - display the processed image (for the 4 quadrants)
;  nocosmic       - cosmic ray/additional bad pixel rejection

function pharoreadfits, fn, header, heap, NOSCALE = noscale, $
  SILENT = silent, EXTEN_NO = exten_no, NUMROW = numrow, $
  POINTLUN = pointlun, STARTROW = startrow, $
  NSLICE = nslice, $
  NO_UNSIGNED = no_unsigned, $
  ;; pharo specialty keywords:
  nocds = nocds, specklecube=specklecube, quadrants=quadrants, $
  ignorenpausefr=ignorenpausefr, refbias=refbias, $
  flatf = flatf, badpixf=badpixf, fixpix_rs=fixpix_rs,  $
  firstbad = firstbad, display = display,  nocosmic = nocosmic

if N_params() LT 1 then begin  
  print,'readfits-pharo readfits wrapper, neads better help'              
  print,'Syntax - im = READFITS( filename, [ h, heap, /NOSCALE, /SILENT,'
  print,'                 NaNValue = ,EXTEN_NO =, STARTROW = , NUMROW='
  print,'                 NSLICE = , /No_UnSigned]'
  print,'                 /NOCDS'
  print,' and more... bug jpl for better docs'
  return, -1
endif

if (keyword_set(firstbad) eq 0) then firstbad =  0
if (keyword_set(nocosmic) eq 0) then nocosmic =  0
if (keyword_set(silent) eq 0) then silent = 0
if (keyword_set(specklecube) eq 0) then specklecube = 0

print,  'About to read: ',  fn
frames=readfits(fn, header, heap, NOSCALE = noscale, $
                SILENT = silent, EXTEN_NO = exten_no, NUMROW = numrow, $
                POINTLUN = pointlun, STARTROW = startrow, $
                ;; NaNvalue = NaNvalue, $
                NSLICE = nslice, $
                NO_UNSIGNED = no_unsigned)

;; ----------------------------------------------------------------------------
;;             substract the dark (if the dark file exists)
;; ----------------------------------------------------------------------------
if keyword_set(refbias) then begin
  info = file_info(refbias)

  ;; ------- refbias may be a directory -------
  if info.directory then $
    refbias = info.name + '/' + pharorefbiasname(header)

  ;; ------- now check whether the file exists or not -------
  info = file_info(refbias)

  if info.exists then begin
    ;; ----- it exists: no problem -----
    frames = frames - readfits(refbias)
    sxaddhist,'refbias used: ' + refbias,header
    print, 'refbias used: ' + refbias
  endif else begin
    ;; ----- it doesn't exist -----
    print, 'WARNING refbias ' + refbias + ' doesnt exist'
    print, 'check parameters or create the dark file (use pharomkrefbias.pro)'
  endelse
endif

writemod=sxpar(header,'WRITEMOD')
naxis1=sxpar(header,'naxis1')
naxis2=sxpar(header,'naxis2')
naxis3=sxpar(header,'naxis3')

;; ----------------------------------------------------------------------------
;; check npausefr first: if != 0, force fowler sampling, whatever specklecube
;; ----------------------------------------------------------------------------

npausefr = sxpar(header, 'npausefr')
if (npausefr ne 0) then begin
  print, 'WARNING NPAUSEFR != 0 : exposures will not be of uniform time'
  print, 'WARNING NPAUSEFR != 0 : Fowler sampling forced'
  specklecube = 0
endif

;; ----------------------------------------------------------------------------
;;   writemod = 0 is the "standard" situation where all frames have been saved
;; ----------------------------------------------------------------------------

if writemod eq 0 then begin

  ;; -----------------------------------------------------------------
  ;;                    pairwise substraction
  ;; -----------------------------------------------------------------
  if keyword_set(specklecube) then begin 
    nframesin = naxis3/4

    ;; -------------------------------------------
    ;;    don't substract cumulative frames
    ;; -------------------------------------------
    if keyword_set(nocds) then begin
      nframesout = nframesin

      ;; -----------------------------------------
      ;;   one keeps the quadrant structure
      ;; -----------------------------------------
      if keyword_set(quadrants) then begin
        frameout = fltarr(naxis1,naxis2, nframesout,4)
        for i=0,nframesin-1 do begin
          for j=0,3 do frameout[*,*,i,j] = frames[*,*,i*4+j]
        endfor                 

      ;; -----------------------------------------
      ;;         quadrants are merged
      ;; -----------------------------------------
      endif else begin
        frameout = fltarr(naxis1*2,naxis2*2, nframesout)
        signal = fltarr(naxis1*2,naxis2*2)
        for i=0,nframesin-1 do begin
          signal = pharojoin(frames[*,*,i*4:i*4+3])
          frameout[*,*,i] = signal
        endfor                 
      endelse

      ;; -------------------------------------------
      ;;       subtract cumulative frames
      ;; -------------------------------------------
    endif else begin
      nframesout = nframesin-1-firstbad

      ;; -----------------------------------------
      ;;   one keeps the quadrant structure
      ;; -----------------------------------------
      if keyword_set(quadrants) then begin
        
        frameout = fltarr(naxis1,naxis2, nframesout,4)
        signal = fltarr(naxis1,naxis2,4)
        reference = fltarr(naxis1,naxis2,4)
        
        for j=0,3 do reference[*,*,j] = frames[*,*,firstbad*4+j]
        
        for i=1+firstbad,nframesin-1 do begin
          for j=0,3 do signal[*,*,j] = frames[*,*,i*4+j]     
          frame = signal - reference
          for j=0,3 do begin
           im =  frame[*,*,j] 
           ;%%%%%  Remove cosmic rays or residual bad pixels (MJI) %%%%%
           if (nocosmic eq 0) then begin
           ; 6-sigma is lots, but remember that we have 65000 pixels.
            im = sigma_filter_nirc2(im,5,n_sigma=5.5,/all,/iterate,  monitor = (silent eq 0))
           ; Remove any 'blinking' pixels not taken into account with bad pixmap
            min_poss = median(im)-3*stdev(im)
            neg =  where(im lt min_poss, nneg)
            if (neg[0] ne -1) then begin
             im[neg] = min_poss
             if (silent eq 0) then print,  nneg, ' residual blinking pixels fixed.'
            endif
           endif      
           frameout[*,*,i-1-firstbad,j] = im     
          endfor
          reference = signal  
          
        endfor     
        ;Display the quadrants.
        if (keyword_set(display)) then begin
         !p.multi = [0, 2, 2]
         image_cont,  total(frameout[*, *, *, 0], 3),  /noc
         image_cont,  total(frameout[*, *, *, 1], 3),  /noc
         image_cont,  total(frameout[*, *, *, 2], 3),  /noc
         image_cont,  total(frameout[*, *, *, 3], 3),  /noc
         !p.multi =  0
       endif
      ;; -----------------------------------------
      ;;         quadrants are merged
      ;; -----------------------------------------
      endif else begin
        frameout = fltarr(naxis1*2,naxis2*2, nframesout)
        reference = fltarr(naxis1*2,naxis2*2)
        signal = fltarr(naxis1*2,naxis2*2)
        
        reference = pharojoin(frames[*,*,0:3])
        
        for i=1,nframesin-1 do begin
          signal = pharojoin(frames[*,*,i*4:i*4+3])
          frameout[*,*,i-1] = signal-reference      
          reference=signal         
        endfor
      endelse
    endelse 

  endif else begin
    ;; --------------------------------------------------------------
    ;;                      Fowler sampling
    ;; --------------------------------------------------------------
    nframesin = naxis3/4
    nsamples  = nframesin/2
    ;;Make sure that firstbad is sensible
    if (firstbad ge nsamples) then begin
        print, 'Firstbad would reject all subframes in the cube! Please change...'
        stop
    endif

    thisframe = fltarr(naxis1, naxis2, 4)
    reference = fltarr(naxis1, naxis2, 4)
    signal    = fltarr(naxis1, naxis2, 4)
    nbS = 0 & nbR = 0 ;; nb of frames in Signal and Ref
    ;; ------------------------------------------------------
    ;; this is a modified version of the Fowler in
    ;; which the frames are weighted (more on extremities).
    ;; ------------------------------------------------------
    for i=nsamples,2*nsamples-1 do begin
      for j=0,3 do thisframe[*,*,j] = frames[*,*,i*4+j]
      signal += (i-nsamples+1) * thisframe
      nbS += i-nsamples+1
    endfor
    signal = signal / float(nbS)

    for i=firstbad, nsamples-1 do begin
      for j=0,3 do thisframe[*,*,j] = frames[*,*,i*4+j]
      reference += (nsamples-i) * thisframe
      nbR += (nsamples-i)
    endfor
    reference = reference / float(nbR)

    ;; ------------------------------------------------------
    ;; now, check for the keywords nocds and quadrants
    ;; ------------------------------------------------------

    ;; -----------------------------------------
    ;;       keep the quadrants structure
    ;; -----------------------------------------
    if keyword_set(quadrants) then begin

      ;; --------------------------------------
      ;;      keep the two average frames
      ;; --------------------------------------
      if keyword_set(nocds) then begin
        nframesout = 2
        frameout = fltarr(naxis1, naxis2, nframesout, 4)
        for j=0,3 do frameout[*,*,0,j]=reference[*,*,j]
        for j=0,3 do frameout[*,*,1,j]=signal[*,*,j]

      ;; --------------------------------------
      ;;   substract the two average frames
      ;; --------------------------------------
      endif else begin
        nframesout = 1
        frameout = fltarr(naxis1, naxis2, 4)
        for j=0,3 do begin
           im = signal[*,*,j]-reference[*,*,j]
           ;%%%%%  Remove cosmic rays or residual bad pixels (MJI) %%%%%
           if (nocosmic eq 0) then begin
           ; 6-sigma is lots, but remember that we have 65000 pixels.
            im = sigma_filter_nirc2(im,5,n_sigma=5.5,/all,/iterate,  monitor = (silent eq 0))
           ; Remove any 'blinking' pixels not taken into account with bad pixmap
            min_poss = median(im)-3*stdev(im)
            neg =  where(im lt min_poss, nneg)
            if (neg[0] ne -1) then begin
             im[neg] = min_poss
             if (silent eq 0) then print,  nneg, ' residual blinking pixels fixed.'
            endif
           endif      
           frameout[*,*,j]= im
        endfor     
      endelse 
       ;Display the quadrants.
        if (keyword_set(display)) then begin
         !p.multi = [0, 2, 2]
         image_cont,  frameout[*, *, 0],  /noc
         image_cont,  frameout[*, *, 1],  /noc
         image_cont,  frameout[*, *, 2],  /noc
         image_cont,  frameout[*, *, 3],  /noc
         !p.multi =  0
        endif
    endif else begin

    ;; -----------------------------------------
    ;;         quadrants are merged
    ;; -----------------------------------------      
      ref = pharojoin(reference)
      sig = pharojoin(signal)
      
      ;; --------------------------------------
      ;;      keep the two average frames
      ;; --------------------------------------
      if keyword_set(nocds) then begin
        nframesout = 2
        frameout = fltarr(2 * naxis1, 2 * naxis2, 2)
        frameout[*,*,0] = ref
        frameout[*,*,1] = sig

      ;; --------------------------------------
      ;;   substract the two average frames
      ;; --------------------------------------
      endif else begin
        nframesout = 1
        frameout = sig - ref
      endelse
    endelse
  endelse
endif else begin

  ;; --------------------------------------------------------------------------
  ;; writemod = 2, "non standard situation": the file only contains 2 sets of
  ;; quadrants frames. It is as if after "Fowler Sampling, nocds, quadrants"
  ;;
  ;; Options here are therefore somewhat limited: merge the quadrants or not !
  ;; --------------------------------------------------------------------------
  print, 'PHARO WRITEMOD = 2'
  
  signal = fltarr(naxis1, naxis2, 4)
  for j = 0, 3 do signal[*,*,j] = frames[*, *, j]
  
  nframesout = 1
  
  ;; --------------------------------------
  ;;     keep the quadrant structure
  ;; i.e.: we do nothing to the data ...
  ;; --------------------------------------
  if keyword_set(quadrants) then begin
    frameout = frames
    
    ;; --------------------------------------
    ;;         merge the quadrants
    ;; --------------------------------------
  endif else begin
    frameout = fltarr(2 * naxis1, 2 * naxis2)
    frameout[*,*] = pharojoin(signal)
  endelse
  
endelse

if nframesout eq 1 then frameout = reform(frameout)

;; ----------------------------------------------------------------------------
;; the pharo detector is a 1024 x 1024: that makes 4 quadrants of 512 x 512
;; now, the fitsframes that this program is reducing are likely to be sub
;; arrays: 256,128 ...
;; ----------------------------------------------------------------------------

;; ----------------------------------------------------------------------------
;;                      data reduction: flat-field
;; ----------------------------------------------------------------------------

if keyword_set(flatf) then begin
  info = file_info(flatf)
  if info.exists then begin
    flat = readfits(flatf)
    ;; check for the bad pixels of the flat field
    dummy = where(flat gt 0.5 and flat lt 2, complement=w)
    flat[w]=1.0  ;; those are already in the badpix map
    flatq = pharosplit(flat)
    
    ;; ---------------- chop the flat down ------------------
    flatq2 = fltarr(naxis1, naxis2, 4)
    for j=0,3 do flatq2[*, *, j] = flatq[0:naxis1-1, 0:naxis2-1, j]
    
    ;; -- we've kept the 4 quadrants structure --
    if keyword_set(quadrants) then begin
      ;; -- how many sets of quadrants do we have ?--
      for i = 0, nframesout - 1 do begin
        if nframesout gt 1 then begin
          for j = 0, 3 do begin
            frameout[*,*,i,j] = frameout[*,*,i,j]/flatq2[*,*,j]
          endfor
        endif else begin
          for j=0,3 do begin
            frameout[*,*,j] = frameout[*,*,j]/flatq2[*,*,j]
          endfor
        endelse
      endfor
      
    endif $
      ;; -- merge the quadrants of the flat then divide --
    else begin
      flat = pharojoin(flatq2)
      for i=0,nframesout-1 do begin
        if nframesout gt 1 then begin
          frameout[*,*,i] = frameout[*,*,i]/flat
        endif
        if nframesout eq 1 then frameout = frameout / flat
      endfor
    endelse
    sxaddhist,'flatf used: ' + flatf, header
    print, 'flatf ' + flatf + ' was used'
  endif else begin
    print, 'WARNING: flat ' + flatf + ' doesnt exist'
    print, 'Data was not flat-fielded'
  endelse
endif

;; ----------------------------------------------------------------------------
;;                    data reduction: bad pixels map
;; ----------------------------------------------------------------------------

if keyword_set(badpixf) then begin
  info = file_info(badpixf)

  if info.exists then begin
    badpix = readfits(badpixf)
    badpixq = pharosplit(badpix)
    
    ;; ----------------- chop the pixmap down ---------------------
    badpixq2 = fltarr(naxis1,naxis2,4)
    for j=0,3 do badpixq2[*,*,j] = badpixq[0:naxis1-1, 0:naxis2-1, j]
    
    if keyword_set(quadrants) then begin
      for i=0,nframesout-1 do begin
        if nframesout gt 1 then begin
          for j=0,3 do begin
            fixpix, frameout[*,*,i,j], badpixq2[*,*,j], temp,  /silent
            frameout[*,*,i,j] = temp
          endfor
        endif else begin
          for j=0,3 do begin
            fixpix, frameout[*,*,j], badpixq2[*,*,j], temp,  /silent
            frameout[*,*,j] = temp
          endfor
        endelse
      endfor
    endif else begin
      ;badpix = pharojoin(badpixq2)
      for j=0,nframesout-1 do begin
        fixpix, frameout[*,*,j], badpix, temp,  /silent
        frameout[*,*,j] = temp
      endfor
    endelse
    sxaddhist,'badpixf used: ' + badpixf,header
  endif else begin
    print, 'WARNING: badpixmap ' + badpixf + ' doesnt exist'
    print, 'Bad pixel treatment has not been performed'
  endelse
endif

;; ----------------------------------------------------------------------------
;;                data reduction: fix bad pixels procedure
;; ----------------------------------------------------------------------------

if keyword_set(fixpix_rs) then begin
  if keyword_set(quadrants) then begin
    for i=0,nframesout-1 do begin
      if nframesout gt 1 then begin
        for j=0,3 do begin
          frameout[*,*,i,j] = fixpix_rs(frameout[*,*,i,j])
        endfor
      endif else begin
        for j=0,3 do begin
          frameout[*,*,i,j] = fixpix_rs(frameout[*,*,i,j])
        endfor
      endelse
    endfor
  endif else begin
    for i=0,nframesout-1 do begin
      frameout = fixpix_rs(frameout)
    endfor
  endelse
  sxaddhist,'fixpix_rs used on all frames'
endif

sxaddhist, 'READFITSPHARO.PRO: '+systime()+userid(), header
sxaddhist, '$Id: pharoreadfits.pro,v 1.11 2009/09/28 23:54:17 snert Exp $', header
sxaddhist, 'fn='+fn, header

;;This is a simple option where read times are doubled, but we end up
;;with only one less frame. Errors remain as independent as
;;specklecube in the case of read-noise dominance only.
if (specklecube eq 2) then begin
 if keyword_set(quadrants) then frameout = frameout[*, *, 1:*, *]+frameout[*, *, 0:nframesout-2, *]$
 else frameout = frameout[*, *, 1:*]+frameout[*, *, 0:nframesout-2]
endif

;;The strange formula below is equivalent to:
;; 2*r5+r4-r2-2*r1 for reads r1, r2, r3, r4, r5. But the formula is
;; written in terms of differences, not reads.
if (specklecube eq 4) then begin
 newnframes = (nframesout-2)/2
 ix = indgen(newnframes)*2-1
 if keyword_set(quadrants) then begin
  frameout = 2*frameout[*, *, 3:*, *]+3*frameout[*, *, 2:nframesout-2, *] + $
   3*frameout[*, *, 1:nframesout-3, *]+ 2*frameout[*, *, 0:nframesout-4, *]
  frameout = frameout[*, *, ix, *]
 endif else begin
  frameout = 2*frameout[*, *, 3:*]+3*frameout[*, *, 2:nframesout-2] + $
   3*frameout[*, *, 1:nframesout-3]+ 2*frameout[*, *, 0:nframesout-4]
  frameout = frameout[*, *, ix] 
 endelse
endif


return, frameout

end

; $Log: pharoreadfits.pro,v $
; Revision 1.11  2009/09/28 23:54:17  snert
; Added features to calibrate_v2_cp suitable for calibrating lots of files.
;
; Revision 1.10  2007/08/02 23:55:05  mireland
; Added specklecube=2 and specklecube=4 modes.
;
; Revision 1.9  2007/06/09 22:31:43  mireland
; Lots of bug fixes and minor feature adds (mainly for using pharomkcube for full
; pupil image analysis).
;
; Revision 1.8  2007/06/02 10:14:00  frantzm
; - Weighted Fowler mode implemented in pharoreadfits.
; - Background substraction (use 3 other quadrants)
;
; Revision 1.7  2007/05/18 17:47:30  mireland
; Fixed a bug with firstbad and added more to the 'silent' option.
;
; Revision 1.6  2006/11/30 16:54:38  frantzm
; just trying to solve a problem with capital letters in the name of a file
; step 1/2
;
; Revision 1.5  2006/09/20 21:53:12  mireland
; Added a display option, a nocosmic and a firstbad option to pharoreadfits.pro.
; Added a specklecube option to the top of the cube_pharo file.
;
; Revision 1.4  2006/02/21 14:32:09  frantzm
; fixed a bug in quadrant flagging: you can now have more cal than src
;
; Revision 1.3  2006/02/09 22:23:24  frantzm
; pharoreadfits now uses a correct bad pixel treatement procedure
;
; Revision 1.2  2006/01/27 16:33:40  frantzm
; check for the existence of auxilliary frames (dark, flat, badpixmap)
;
; Revision 1.1.1.1  2005/12/19 05:15:05  mireland
; Initial release of masking code.
;
; Revision 1.2  2003/12/31 02:20:02  jpl
; continued development, including flatfields and bad pixel treatment
; for pharoreadfits
;
; Revision 1.1  2003/12/19 23:03:40  jpl
; *** empty log message ***
;
; Revision 1.8  2003/10/28 02:10:29  jpl
; tweaks and twiddles in Sydney, Oct 03
;
; Revision 1.7  2003/10/24 07:35:18  jpl
; change to be backward compatible with idl_5.6
;
; Revision 1.6  2003/10/24 06:47:21  jpl
; fix fowler sampling
;
; Revision 1.5  2003/10/24 02:16:19  jpl
; *** empty log message ***
;
; Revision 1.3  2003/09/16 20:46:27  jpl
; fix the fowler sampling - largely untested
;
; Revision 1.2  2003/09/14 11:02:43  jpl
; add specklecube mode

