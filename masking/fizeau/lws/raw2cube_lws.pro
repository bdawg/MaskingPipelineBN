;+
; NAME:
;      raw2cube_lws.pro
; PURPOSE:
; 	Analyze multiple (rapid) LWS frames 
; 	Does sky subtraction/bad pixel correction/flat fielding. 
;
; USAGE:
;   raw2cube_lws,
;
; INPUTS:
;	frame_num - file number of frame to be analyzed
;
; OPTIONAL INPUTS:
;	gain	- gain of chip (flat field)
;	badpix  - vector of bad pixels
;	datadir - path where data is kept (default ./)
;       showraw - =0  default operation
;                 =1  show each frame as it is shifted!
;                 =2  Click to accept/reject each speckle frame
;                 =3  Click on centroid of each speckle frame
;
; OUTPUTS:
;cleancube	- output cleaned frames
;    ROTN	- rotator info from header (Structure)
;    eADU	- electron/ADU from header 
;
; Written					28Aug04 PGT
;  frame skipping (reject) code does not work yet...

pro raw2cube_lws,frame_num,cleancube,nodsky=nodsky, $
          gain=gain,badpix=badpix,datadir=datadir, $
          ROTN=ROTN,eADU=eADU,showraw=showraw,  nod_ims = nod_ims,  head = head,  fix_weird_fits = fix_weird_fits

if (keyword_set(frame_num) eq 0) then begin
 print,'USAGE:'
 print,'   raw2cube_lws.pro,frame_num,cleancube'
 return
endif

if (keyword_set(datadir) eq 0) then datadir='./'

if (keyword_set(showraw) eq 0) then showraw=0  ; 0=regular operation
if(showraw eq 2) then begin 
  skipcount=0
  print,''
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,'You have selected the option to manually discard frames'
  print,'To accept a frame, click right of the y axis'
  print,'To reject a frame, click left  of the y axis'
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,''
endif
if(showraw eq 3) then begin
  skipcount=0
  print,''
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,'You have selected the option to manually Center frames'
  print,'Click on the centroid of each frame'
  print,'To reject a frame, click left  of the y axis'
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,''
endif


; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; First part of code - deal with Nod sky frames 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (keyword_set(nodsky) ne 0) then begin

; Read in data
  filename=datadir+"lws"+string(nodsky,format="(I4.4)")+".fits"
  c=float(readfits(filename,head))
  sz=size(c)

NAX  = sxpar(head,'NAXIS')    ; number of axes            
AX1  = sxpar(head,'NAXIS1')   ; COLUMNS (128) 
AX2  = sxpar(head,'NAXIS2')   ; ROWS    (128)   
AX3  = sxpar(head,'NAXIS3')   ; CHOP BEAMS 
AX4  = sxpar(head,'NAXIS4')   ; CHOP SETS  
AX5  = sxpar(head,'NAXIS5')   ; NOD BEAMS  
AX6  = sxpar(head,'NAXIS6')   ; FINAL NOD SETS
frmcoadd=sxpar(head,'FRMCOADD') 
chpcoadd=sxpar(head,'CHPCOADD')
ncoadd=chpcoadd*frmcoadd


; If we used the fastspk timing pattern, the data are inverted
; and we need to subtract sky-source ... use the variable ADU_gain
; in order to achieve this 
  eADU=float(strtrim(sxpar(head,'GAIN'),2))   ; electrons/ADU
  if(eADU eq 600.0) then ADU_gain=-1.0 else ADU_gain=1.0

  if (keyword_set(gain) eq 0) then begin
   gain=replicate(1.0,AX1,AX2)
  endif

  if(sz(0) lt 4) then begin
    print,'### Unexpected number of dimensions in file ###'
    stop
  endif
  if(AX3 lt 2) then begin
    print,'### WOT no Chopping ?? ###'
    stop
  endif

; Firstly sum over chops sets:
  odds =total(c(*,*,0,*,*,*),4)
  evens=total(c(*,*,1,*,*,*),4)

; Now (if necessary) sum over nod beams/sets
  if(AX5 eq 2) then begin       ;we have 2 nod beams
    print,'### This is pretty bizzare - I was not expecting multiple Nods! ###'
    if(AX6 gt 1) then begin     ;we have >1 nod set
       fsrc=total(odds(*,*,0,1,*),5)+total(evens(*,*,0,0,*),5)
       fsky=total(odds(*,*,0,0,*),5)+total(evens(*,*,0,1,*),5)
    endif else begin	      ;only 1 nod set
       fsrc=odds(*,*,0,1)+evens(*,*,0,0)
       fsky=odds(*,*,0,0)+evens(*,*,0,1)
    endelse
  endif else begin              ; NO nodding at all
    fsrc=reform(evens)
    fsky=reform(odds)
  endelse
  
  nodframe=reform(fsrc-fsky)/float(ncoadd * AX4 * AX5 * AX6) * ADU_gain
  nodframe=nodframe/gain
  if (keyword_set(badpix)) then $
   nodframe=fix_bad_pixels(nodframe,bad=badpix)
 ;--- Mike's new test stuff ---
 nod_ims =  fltarr(sz[1], sz[2],  sz[4])
 for i = 0,  sz[4]-1 do $
  nod_ims[*, *, i]= reform(c[*,*,1,i] - c[*,*,0,i])/gain/float(ncoadd * AX5 * AX6)* ADU_gain - nodframe
 if (keyword_set(badpix)) then for i = 0,  sz[4]-1 do $
   nod_ims[*, *, i] = fix_bad_pixels(nod_ims[*, *, i], bad = badpix)
endif else nodframe=-1

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Now read in and process Data frames
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; Read in data
filename=datadir+"lws"+string(frame_num,format="(I4.4)")+".fits"
;c=float(reform(readfits(filename,head)))
c=readfits(filename,head)
if (keyword_set(fix_weird_fits)) then c =  uint(c) ;!!!Unsure if this will always work...
c = float(c)

NAX  = sxpar(head,'NAXIS')    ; number of axes            
AX1  = sxpar(head,'NAXIS1')   ; COLUMNS (128) 
AX2  = sxpar(head,'NAXIS2')   ; ROWS    (128)   
AX3  = sxpar(head,'NAXIS3')   ; CHOP BEAMS 
AX4  = sxpar(head,'NAXIS4')   ; CHOP SETS  
AX5  = sxpar(head,'NAXIS5')   ; NOD BEAMS  
AX6  = sxpar(head,'NAXIS6')   ; FINAL NOD SETS

if (n_elements(nodframe) eq 1) then nodframe =  fltarr(AX1, AX2)
cleancube=fltarr(ax1,ax2,ax4*ax5)
frmcoadd=sxpar(head,'FRMCOADD') 
chpcoadd=sxpar(head,'CHPCOADD')
; chpsets =sxpar(head,'SAVESETS')
; chpbeams=sxpar(head,'CHPBEAMS')
; nodbeams=sxpar(head,'NODBEAMS') 
; nodsets =sxpar(head,'NODSETS') 
ncoadd=chpcoadd*frmcoadd

rotn={PARANTEL:0.0,ROTPOSN:0.0,ROTPPOSN:0.0,ROTMODE:''}
rotn.PARANTEL=sxpar(head,'PARANTEL')
rotn.ROTPOSN=sxpar(head,'ROTPOSN')
rotn.ROTPPOSN=sxpar(head,'ROTPPOSN')
rotn.ROTMODE=sxpar(head,'ROTMODE')

; If we used the fastspk timing pattern, the data are inverted
; and we need to subtract sky-source ... use the variable ADU_gain
; in order to achieve this 
eADU=float(strtrim(sxpar(head,'GAIN'),2))   ; electrons/ADU
if(eADU eq 600.0) then ADU_gain=-1.0 else ADU_gain=1.0

if (keyword_set(gain) eq 0) then begin
 gain=replicate(1.0,AX1,AX2)
endif

sz=size(c)

; Caution! IDL will discard dimensions of 1 BUT ONLY AT THE END
; thus array(2,3,2,1,1,) --> array(2,3,2) BUT array(2,1,2) --> array(2,1,2)

if(sz(0) lt 4) then begin
  print,'### Unexpected number of dimensions in file ###'
  stop
endif
if(AX3 lt 2) then begin
  print,'### WOT no Chopping ?? ###'
  stop
endif


;CASE1: Multiple chop sets, Nodding OFF
if(AX4 gt 1 and AX5 eq 1 and AX6 eq 1) then begin 
print,AX4,' Chop Sets. No Nodding'
   for cs=0,AX4-1 do begin
         spk=(c(*,*,1,cs)-c(*,*,0,cs)) / ncoadd * ADU_gain
         spk=spk/gain-nodframe ;nodframe has already been divided by gain...
         spk=fix_bad_pixels(spk,bad=badpix)
         ;  if(showraw gt 0) then begin  ;Manual Frame operation
         ;    image_cont,spk,/nocont,/asp 
         ;    if(showraw ge 2) then begin 
         ;      cursor,x,y
         ;      if(x lt 0) then begin
         ;        skipcount=skipcount+1
         ;        goto,skipframe1
         ;      endif 
         ;      if(showraw eq 3) then begin & x0=x & y0=y & endif
         ;    endif else wait,.5
         ;  endif                   ; END Manual Frame operation
         ;  skipframe1:
         cleancube[*,*,cs]=spk
   endfor
endif else if(AX4 gt 1 and AX5 gt 1 and AX6 eq 1) then begin 
;CASE2: Multiple chop sets, Nodding ON, one NOD set only
print,AX4,' Chop Sets. One Nod Pair'
for nod=0,AX5-1 do begin
   for cs=0,AX4-1 do begin
      if(2*(nod/2) eq nod) then begin
         spk=(c(*,*,1,cs,nod)-c(*,*,0,cs,nod)) / ncoadd * ADU_gain
         spk=spk/gain-nodframe
         spk=fix_bad_pixels(spk,bad=badpix)
      endif else begin
         spk=(c(*,*,0,cs,nod)-c(*,*,1,cs,nod)) / ncoadd * ADU_gain
         spk=spk/gain+nodframe
         spk=fix_bad_pixels(spk,bad=badpix)
      endelse
      cleancube[*,*,cs+nod*AX4]=spk
   endfor
endfor
endif else print,'CHOP-NOD sets in this file not dealt with by this program'

if(showraw ge 2) then begin 
  print,''
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,'You have rejected',skipcount,' frames'
  print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print,''
endif

;stop

end

