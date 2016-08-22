;+
; NAME:
;      raw2cube_trecs.pro
; PURPOSE:
; 	Analyze multiple (rapid) TReCS frames 
; 	Does nod subtraction/bad pixel correction/trimming off array
;
; USAGE:
;   raw2cube_trecs,
;
; INPUTS:
;	frame_num - file number of frame to be analyzed
;
;!!! Stuff below here to be changed !!!
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

pro raw2cube_trecs, frame_num,cleanframes,nodsky=nodsky,  bad_pixels = bad_pixels,  destripe = destripe, $
          gain=gain,datadir=datadir, findbad = findbad, $
          ROTN=ROTN,eADU=eADU,showraw=showraw,  nod_ims = nod_ims,  head = head

;if (keyword_set(cleancube) eq 0) then begin
; print,'USAGE:'
; print,'   raw2cube_trecs.pro, frame_num,cleancube'
; return
;endif

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



; Read in data. NB datadir includes the file prefix
  filename=datadir+string(frame_num,format="(I4.4)")+".fits.gz"
  dummy=float(readfits(filename,mainhead))
  nnods = sxpar(mainhead, 'NNODSETS') ;Nuber of nods. Twice this number of extensions
  dummy =  reform(float(readfits(filename,head,  exten = 1)))
 sz=size(c)

NAX  = sxpar(head,'NAXIS')    ; number of axes            
AX1  = sxpar(head,'NAXIS1')   ; COLUMNS (320) 
AX2  = sxpar(head,'NAXIS2')   ; ROWS    (240)   
AX3  = sxpar(head,'NAXIS3')   ; ??? 
AX4  = sxpar(head,'NAXIS4')   ; Number of frames

 c =  fltarr(ax1, ax2, ax4*nnods, 2)
 for i = 1, nnods*2 do begin
   if (i mod 2) eq 1 then $
    c[*, *, (i-1)/2*ax4:(i+1)/2*ax4-1, (i-1) mod 2]=  reform(float(readfits(filename,head,  exten = i))) $
   else $
    c[*, *, (i-1)/2*ax4:(i+1)/2*ax4-1, (i-1) mod 2]=  reform(float(readfits(filename,head,  exten = i)))
 endfor

eADU=100.0;  !!! Just a guess !!! float(strtrim(sxpar(head,'GAIN'),2))   ; electrons/ADU

;Bad Pix
if (keyword_set(bad_pixels)) then for j = 0, 1 do for i = 0, ax4*nnods-1 do begin
 c[*, *, i, j] =  fix_bad_pixels(c[*, *, i, j], bad_pix = bad_pixels)
endfor
;De-striping
if (keyword_set(destripe)) then for j = 0, 1 do for i = 0, ax4*nnods-1 do begin
 line = total(c[*,10:29,i,j]+c[*,220:239,i,j],2)/40.
 for k = 0, 239 do c[*, k, i, j] -= line
; image_cont,  c[*, *, i, j],  /noc &  wait,  0.2
endfor

if (keyword_set(findbad)) then begin
 medim = total(median(c, dim=3), 3)/2.0
 mn = fltarr(ax1, ax2)
 mn2 =  fltarr(ax1, ax2)
 for j = 0, 1 do for i=0,ax4*nnods-1 do mn2 += (c[*, *, i, j]-medim)^2
 mn2 /= 2*ax4*nnods
 mn2 /= mean(mn2)
 newbad =  where(mn2 gt 5)
 newbad =  [bad_pixels, newbad]
 newbad =  newbad[uniq(newbad)]
 print,  'New Bad pixels stored in newbad. You should save them...'
 stop
endif

mnim = total(c, 3)/ax4/nnods
subim =  mnim[*, *, 0]-mnim[*, *, 1]
pos =  fltarr(2, 2)
dummy = max(subim,m)
pos[*, 0] =  array_indices(subim, m)
dummy = min(subim,m)
pos[*, 1] =  array_indices(subim, m)

cleanframes =  fltarr(128, 128, ax4*nnods*2)
d = shift(dist(128), 64, 64)
annpix =  where(d gt 50 and d lt 70)
for i = 0, ax4*nnods-1 do begin
; cleanframes[*, *, i] = c[pos[0, 0]-64:pos[0, 0]+63, pos[1, 0]-64:pos[1, 0]+63, i, 0] - $
;   mnim[pos[0, 0]-64:pos[0, 0]+63, pos[1, 0]-64:pos[1, 0]+63, 1]
; cleanframes[*, *, i+ax4*nnods] = c[pos[0, 1]-64:pos[0, 1]+63, pos[1, 1]-64:pos[1, 1]+63, i, 1] - $
;   mnim[pos[0, 1]-64:pos[0, 1]+63, pos[1, 1]-64:pos[1, 1]+63, 0] 
 temp =  (shift(c[*, *, i, 0]-mnim[*, *, 1], -pos[0, 0]+64, -pos[1, 0]+64))[0:127, 0:127]
 cleanframes[*, *, i] = temp - median(temp[annpix])
 temp = (shift(c[*, *, i, 1]-mnim[*, *, 0], -pos[0, 1]+64, -pos[1, 1]+64))[0:127, 0:127]
 cleanframes[*, *, i+ax4*nnods] = temp-median(temp[annpix])
endfor
image_cont,  subim,  /noc;total(cleanframes, 3),  /noc
head =  mainhead

end

