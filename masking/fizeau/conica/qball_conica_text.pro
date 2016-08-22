;######## RELEASE CLIP HERE ######## 
;+
;qball.script
;
; A tool to assist making first pass cube-all data in an entire directory.
; Pulls useful information from all headers of a set of data files. 
; and tries to produce a set of cubes clumping files together when it
; finds something changes 
; PGT   May 2009



data_dir='~snert/conica/data/conica_Mar10/'
;data_dir='/import/quoll2/snert/conica/data/conica_Mar09/082.C-0742A/2009-03-06/SCIENCE/'

flats_dir= '~snert/conica/data/conica_Mar10/'
;flats_dir='/import/quoll2/snert/conica/data/conica_Mar09/Testing/'

src_names=['SR21']
nsrc_names=n_elements(src_names)
root_dir='~/code/masking/'
prefix='NACO_IMG_SCI076_'
extn='.fits.gz'
headtxt_output='header.txt'


;flags 
dsky           = -1     ; default skies value, best to leave it set to dithering (-1)
                        ;     here and change later with option 3
updown         = 0    	; 0= normal operation no wollaston, find star anywhere on chip
			; 1= find star in top half of chip
			; 2= find star in bottom half of chip
horiz_destripe = 1      ; 1 = basic horizontal de-striping is applied to data after 
                        ;     sky subtraction.
                        ; 0 = no horizontal de-striping applied.                                  
discard_last   = 1      ; 1 = last frame of each cube is not included in the analysis
                        ;     (eg. when last frame is average of all other frames). 
                        ; 0 = last frame included in the analysis. 
displayall     = 0      ; 0=display only changes; 1=print every single file


; %%%%%%%%%%%%%%%%%%%%
; AUTOMATIC BELOW HERE
; %%%%%%%%%%%%%%%%%%%%

; Find the current working directory
CD, CURRENT=original_dir
   starname='starname'
   dt= systime()
   datename=strmid(dt,8,2)+strmid(dt,4,3)+strmid(dt,22,2)

; Housekeeping to save/restore previous qball markups ...
qb_name=''
read, qb_name, prompt='Name of this data markup (e.g. qb_conica_14Jun09.idlvar) > ' 
print,'Restore previous y/N? > '
d = get_kbrd(1) & wait,.5 
;if the file list is not in memory then it needs to be re-read in, as this is not stored in the idlvar above
if(d eq 'y' or d eq 'Y') then if(n_elements(file_id) ne 0) then goto,skipheadread else goto, readhead

readhead:
; Things to grab from header ...
s_par=['ORIGFILE','OOBSNAME','DATE-OBS','SOPTI4ID','SOPTI5ID','SOPTI6ID','SOPTI3ID',  'RETA2ROT', 'SOPTI7ID']  
n_par=['RA','DEC','NAXIS1','NAXIS3','SODETDIT','SOTELALT ','ARANGEND']
 
if (keyword_set(data_dir) eq 0) then data_dir="./"
file=findfile(data_dir+prefix+'*'+extn)
if file[0] eq "" then begin
    print, "No files found"
    goto, last
endif
n_files=n_elements(file)
spars=n_elements(s_par)
npars=n_elements(n_par)
s_output=strarr(spars,n_files)
n_output=fltarr(npars,n_files)
pos_ang=fltarr(n_files)
pos_ang_lessPA=fltarr(n_files)


; Read through all the headers, grab relevant info...
new_src=[0]
print, strcompress(string(n_files)), " files found. Reading..."
for i=0,n_files-1 do begin
   filename=file(i)
   head=headfits(filename)
  for j=0,n_elements(s_par)-1 do begin
   s_output(j,i)=sxpar_conica(head,s_par(j))
  endfor
  for j=0,n_elements(n_par)-1 do begin
   n_output(j,i)=sxpar_conica(head,n_par(j))
  endfor
  ; Formula for position angle taken from /import/spiral1/snert/code/masking/fizeau/freud.pro
  rotstart = sxpar_conica(head,'ROTSTART')
  rotend   = sxpar_conica(head,'BSROTEND')
  pastart  = sxpar_conica(head,'ANGSTART')
  paend    = sxpar_conica(head,'ARANGEND')
  alt      = sxpar_conica(head,'SOTELALT')
  instrument_offset = -0.55
  pos_ang(i) = (rotstart+rotend)/2.+alt-(180.-(pastart+paend)/2.) + instrument_offset
  pos_ang_lessPA(i) = (rotstart+rotend)/2.+alt + instrument_offset 
  ; work out if the present observation is a new target/color/nod or simply another data file same as before ...
    s_ix=[1,3,4,5,6,7,8]  ; string quantities to check
  if(i gt 0) then  begin
    if( abs(n_output(0,i)-run_ra) gt 1. or       $   ; change in pointing?
        abs(n_output(1,i)-run_dec) gt .01 or     $   ;
        strjoin(s_output(s_ix,i)) ne run_strings or $   ; change in config (strings)?
        strjoin(string(n_output(2:4,i))) ne run_ints ) $ ; change in config (int)?
      then new_src = [new_src,i]
  endif
  run_ra=n_output(0,i)
  run_dec=n_output(1,i)
  run_strings=strjoin(s_output(s_ix,i))
  run_ints=strjoin(string(n_output(2:4,i)))
endfor   

if(displayall eq 1) then new_src=indgen(n_files)  ; for now, work with every file individually (don't try to group similar ones)

; Now do pretty tabular output of numbers
; Trim the string array
s_output=strcompress(strtrim(s_output,2))
;s_output(0,*)=strcompress(s_output(0,*),/remove_all)  ; This gets rid of whitespace in filter
s_output=strcompress(s_output,/remove_all)  ; This gets rid of _ALL_ whitespace

n_lines=n_elements(new_src)
s_output=s_output(*,new_src)
n_output=n_output(*,new_src)
filen=file(new_src)
filenumlimits=intarr(2,n_lines)
print,' IX     File#ID           ObsBlk      TIME       RA           DEC    Filter  OPTI4ID Mask    AX1  AX3  T_int  ALT    PosAng PupilRot  HWPlt CAM'
;s_output(0,*) = file ; only needed when there is an error in the embedded filenames
for i=0,n_lines-1 do begin
  pos=strpos(s_output(0,i),strmid(extn,0,5))
  file_id=strmid(s_output(0,i),pos-15,15)
  file_num=fix(strmid(s_output(0,i),pos-4,4))
  filenumlimits[0,i]=file_num
  if(where(new_src eq new_src(i)+1) eq -1 and i lt n_lines-2) then begin
    pos2=strpos(file(new_src(i+1)-1),extn)
    file_interval='-'+ strmid(file(new_src(i+1)-1),pos2-2,2)
    filenumlimits[1,i]=fix(strmid(file(new_src(i+1)-1),pos2-4,4))
  endif else begin
    file_interval=''
    filenumlimits[1,i]=filenumlimits[0,i]
    if(i lt n_lines-2) then begin
       pos2=strpos(file(n_files-1),extn) 
       filenumlimits[1,i]=fix(strmid(file(fix(strmid(s_output(0,i+1),pos-4,4))-1),pos2-4,4))
       file_interval='-'+ strmid(file(n_files-1),pos2-2,2)
    endif
  endelse
  if(i eq n_lines-1 and new_src(i) ne n_files-1) then file_interval='-'+ strmid(file(n_files-1),pos-2,2)
  ; Trim the timestamp information:
  redtime=strmid(s_output(2,i),11,8)
  ; *** find the filter 
  filter = s_output(4,i)  ; i.e. SOPTI5ID
  if (filter eq 'empty') then filter = s_output(5,i)  ; i.e. SOPTI6ID
  if (filter eq 'empty') then filter = s_output(3,i)  ; i.e. SOPTI4ID
  print,i,file_id,s_output(1,i),redtime,adstring(n_output[0,i],n_output[1,i]),filter,s_output(3,i),s_output(6,i), $
         n_output(2:5,i),pos_ang[i],pos_ang_lessPA[i],(360.0+s_output(7,i)) mod 360,s_output(8,i), $
         format="(I3,' ',A15,' ',A15,' ',A8,' ',A22,' ',A8,' ',A6,' ',A7,' ',I4,' ',I4,' ',F6.2,' ',F6.2,' ',F7.2,' ',F7.2,' ',F6.2,' ',A4)"
endfor
if d eq 'y' or d eq 'Y' then goto, skipheadread
; Now write to output file
openw,2,headtxt_output
printf,2,' IX     File#ID           ObsBlk      TIME       RA           DEC    Filter  OPTI4ID Mask    AX1  AX3  T_int  ALT    PosAng PupilRot  HWPlt CAM'
for i=0,n_lines-1 do begin
  pos=strpos(s_output(0,i),'.fits')
  file_id=strmid(s_output(0,i),pos-15,15)
  ; Trim the timestamp information:
  redtime=strmid(s_output(2,i),11,8)
  ; *** find the filter 
  filter = s_output(4,i)  ; i.e. SOPTI5ID
  if (filter eq 'empty') then filter = s_output(5,i)  ; i.e. SOPTI6ID
  if (filter eq 'empty') then filter = s_output(3,i)  ; i.e. SOPTI4ID
  printf,2,i,file_id,s_output(1,i),redtime,adstring(n_output[0,i],n_output[1,i]),filter,s_output(3,i),s_output(6,i), $
         n_output(2:5,i),pos_ang[i],pos_ang_lessPA[i],(360.0+s_output(7,i)) mod 360,s_output(8,i), $
         format="(I3,' ',A15,' ',A15,' ',A8,' ',A22,' ',A8,' ',A6,' ',A7,' ',I4,' ',I4,' ',F6.2,' ',F6.2,' ',F7.2,' ',F7.2,' ',F6.2,' ',A4)"
endfor
close,2

n_clumps=0

skipheadread:
if d eq 'y' or d eq 'Y' then restore,qb_name

selectmode:
print,'  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print,'  % Select option:                                  %'
print,'  % 0 = Print original header strip                 %'
print,'  % 1 = Define a new data Block                     %'
print,'  % 2 = Show presently defined Blocks               %'
print,'  % 3 = Edit a Block                                %'
;Option 4 not working properly, delete then define new works properly,
; but will change order. I don't think that is a problem...
;The problem is when you replace a block with one that has a different size,
; the rest of the blocks get disrupted, either pushed or pulled into other blocks.
;print,'  % 4 = Overwrite a Block                           %'
print,'  % 5 = Delete a Block                              %'
print,'  % 8 = Save present state of data markups          %'
print,'  % 9 = Finished with data tagging - queue analysis %'
print,'  % q = quit                                        %'
print,'  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
k = get_kbrd(1) & wait,.5
case k of


; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 0 %%%%%%%%%%%%%%%%%%%%%%%%%
  '0': begin   ; Print stuff out...
  if file_test(headtxt_output) eq 1 then spawn, 'more '+headtxt_output else begin
print,' IX  Block  File#ID           ObsBlk      TIME       RA           DEC    Filter  OPTI4ID Mask    AX1  AX3  T_int  ALT    PosAng PupilRot  HWPlt CAM'
for i=0,n_lines-1 do begin
  pos=strpos(filen(i),extn)
  file_num=fix(strmid(filen(i),pos-4,4))
  if(where(new_src eq new_src(i)+1) eq -1 and i lt n_lines-1) then begin
    pos2=strpos(file(new_src(i+1)-1),extn)
    file_interval='-'+ strmid(file(new_src(i+1)-1),pos-2,2)
  endif else file_interval=''  
  blockstring=''
  bw=where(flagblock[1,*] eq i) 
    if(bw[0] ne -1) then begin
      if(n_elements(bw) eq 1) then blockstring=strtrim(string(flagblock[2,bw]),2) $
        else print,'* NOTE - this star set in mutiple blocks!'
      bw=bw[0]
    endif
  ;if(i eq n_lines-1 and new_src(i) ne n_files-1) then file_interval='-'+ strmid(file(n_files-1),pos-2,2)
  print,i,blockstring,file_id,s_output(1,i),redtime,adstring(n_output[0,i],n_output[1,i]),filter,s_output(3,i),s_output(6,i), $
         n_output(2:5,i),pos_ang[i],pos_ang_lessPA[i],(360.0+s_output(7,i)) mod 360,s_output(8,i), $
         format="(I3,' ',A4,' ',A15,' ',A15,' ',A8,' ',A22,' ',A8,' ',A6,' ',A7,' ',I4,' ',I4,' ',F6.2,' ',F6.2,' ',F7.2,' ',F7.2,' ',F6.2,' ',A4)"
  ;print,i,blockstring,file_num,file_interval,s_output(*,i),n_output(0:7,i),format="(I4,A4,I4,A3,' ',A17,2(A12),A7,A7,A9,A14,A12,3(I5),I2,F8.3,F5.1,2F8.2)"
endfor  
  endelse
end


; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 1 %%%%%%%%%%%%%%%%%%%%%%%%%
  '1': begin   ; define a new group of files to calibrate together

   starnameo = starname
   datenameo = datename
   starname = ''
   datename = ''
   read, starname, prompt='Star identifier/working dir for the new grouping (default: '+starnameo+') > ' 
   read, datename, prompt='Date identifier for the new grouping (default: '+datenameo+') > ' 
   if starname eq '' then starname = starnameo
   if datename eq '' then datename = datenameo
   repeat begin
   print,'Type the indexes of the stars in star grouping number ',strtrim(string(n_clumps),2), " e.g. 9 10 11 )
   instr=''
   read, instr, prompt='> ' ; required to take multiple integers as input
   thisclump=fix(string_vector(instr))  
   endrep until min(thisclump) ge 0 and max(thisclump) lt n_lines
   n_stars_thisclump=n_elements(thisclump)
   if(n_clumps eq 0) then n_stars_per_clump=[n_stars_thisclump] $
     else n_stars_per_clump=[n_stars_per_clump,n_stars_thisclump]
   for s=0, n_stars_thisclump-1 do begin 
     zvect=intarr(100)
     flagvect=intarr(10)
     nstarfiles=filenumlimits[1,thisclump[s]]-filenumlimits[0,thisclump[s]]+1
     thisfilevect=indgen(nstarfiles) + filenumlimits[0,thisclump[s]]
     zvect[0:nstarfiles-1]=thisfilevect  ; use array elements 10+ for the filenums
      ; Array elements 0-9 are reserved for information...
     flagvect[0]=nstarfiles    ; num files this star
     flagvect[1]=thisclump[s]  ; #index number through night
     flagvect[2]=n_clumps      ; The Block number of this clump
     flagvect[3]=1             ; Default setting source/cal flag
     ; Try to figure out if it is a source or cal
     for t=0,nsrc_names-1 do $
         if( strpos(s_output[1,thisclump[s]],src_names[t]) ne -1) then flagvect[3]=0  ; Prolly a source
     if( strpos(s_output[1,thisclump[s]],'SCI') ne -1) then flagvect[3]=0             ; Prolly a source	 
     if( strpos(s_output[1,thisclump[s]],'CAL') ne -1) then flagvect[3]=1             ; Prolly a CAL	 
     if( strpos(s_output[1,thisclump[s]],'PSF') ne -1) then flagvect[3]=1             ; Prolly a CAL
     flagvect[4]=1             ; 0=nodither; 1=dither
     flagvect[5]=1             ; Observation Type:  1=normal src/cal;   2=faint companion hunt
     ;dithering is on by default, changed in option 3 where skies can be added
     sky = create_struct('path',ptr_new(""),'files', ptr_new(dsky))
     if(s eq 0 and n_clumps eq 0) then begin
          fileblock=zvect 
          flagblock=flagvect
          s_block  =s_output[*,thisclump[s]]
          n_block  =n_output[*,thisclump[s]]
          datenames = datename
          starnames = starname
	  skies = sky
	  bispectflags = ""
	  calv2cpflags = ""
	  bingridflags = ""
	  
     endif else begin
          fileblock=[[fileblock],[zvect]] 
          flagblock=[[flagblock],[flagvect]]
          s_block  =[[s_block],[s_output[*,thisclump[s]]]]
          n_block  =[[n_block],[n_output[*,thisclump[s]]]]
          datenames = [datenames,datename]
          starnames = [starnames,starname]
	  skies = [skies,sky]
	  bispectflags = [bispectflags,""]
	  calv2cpflags = [calv2cpflags,""]
	  bingridflags = [bingridflags,""]
     endelse
   endfor
   n_clumps=n_clumps+1
   end


; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 2 %%%%%%%%%%%%%%%%%%%%%%%%%
  '2': begin   ; show present blocks
   if n_clumps eq 0 then print, "No Currently defined blocks" else begin
   print,'[ Flags Note: SRCCAL 0/1 = Src/Cal || DITHER 0/1 = No/Yes || TYPE 1/2 = normal/faint_companion ]'
   print,''
   print,'|  REF NUMBERS   |       |       FLAGS   Obs  |'
   print,'|Index Block Star N_Files SrcCal Dither  Type | Data Filenumbers...'
   ixcount=0
   for b=0,n_clumps-1 do begin 
     print,'                                                Star Block ',starnames[ixcount],'_',datenames[ixcount]
     for s=0,n_stars_per_clump[b]-1 do begin 
       w=where(fileblock[*,ixcount] gt 0)
       print,flagblock[1,ixcount],flagblock[2,ixcount],s,flagblock[0,ixcount],flagblock[3,ixcount], $
              flagblock[4,ixcount],flagblock[5,ixcount],fileblock[w,ixcount],format="(3I5,5I7,100I5)"
       ixcount += 1 
     endfor
   endfor
   endelse
   end

; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 3 %%%%%%%%%%%%%%%%%%%%%%%%%
  '3': begin   ; Edit blocks
   repeat begin
   print,'Input the Block# of the block to be edited, 0 -'+strcompress(string(n_clumps-1))
   read, instr, prompt='> ' 
   endrep until fix(instr) lt n_clumps and fix(instr) ge 0
   
   thisblock=fix(instr) 
   bb=where(flagblock[2,*] eq thisblock)
   
; modify the DITHER and TYPE flags if requested  ...
   print,'This block: Dither, Type = ',strtrim(string(flagblock[4:5 ,bb[0]]),2),' : modify y/N? > '
   d = get_kbrd(1) & wait,.5 

   if(d eq 'y' or d eq 'Y') then begin
     repeat begin
     print,'Input new values for Dither and Type e.g. > 0  or > 1 2 (Dither 0/1 = Off/On  TYPE 1/2 = normal/faint_companion) ' 
     read, instr, prompt='>' 
     newline=fix(string_vector(instr))   
     endcond = 0
     if n_elements(newline) eq 1 and newline[0] eq 0 then endcond = 1 else $
     if n_elements(newline) eq 2 and newline[0] eq 1 and (newline[1] eq 1 or newline[1] eq 2) then endcond = 1
     endrep until endcond eq 1
       flagblock[4,bb]=newline[0]
       if newline[0] eq 1 then flagblock[5,bb]=newline[1] else flagblock[5,bb]=1
     ;determine the sky files
       if newline[0] eq 0 then begin
           print,'Input path of the sky files, blank for the same as data >' 
           read, instr, prompt='>'        
           *skies[bb[0]].path = instr
           print,'Input prefix of the sky files >' 
           read, instr, prompt='>' 
	   skyprefix = instr
           print,'Input suffix of the sky files >' 
           read, instr, prompt='>' 	   
	   skysuffix = instr
           print,'Input starting file number >' 
           read, instr, prompt='>' 	   
	   skystart = fix(instr)
           print,'Input number of sky files >' 
           read, instr, prompt='>' 	   
	   nskies = fix(instr)
	   
           *skies[bb[0]].files= skyprefix+string(skystart+indgen(nskies),format="(I4.4)")+skysuffix 
       endif
   endif
   print,''
; modify the SRC/CAL flag if requested  ...
   print,'This block Src(0)/Cal(1) = ',strtrim(string(reform(flagblock[3,bb])),2)
   print,'Modify y/N? > '
   d = get_kbrd(1) & wait,.5 
   if(d eq 'y' or d eq 'Y') then begin
     tryscagain:
     print,'Input new values for Src/Cal flag. 0=Src, 1=Cal eg, >0 0 1 1' 
     read, instr, prompt='>' 
     newline=fix(string_vector(instr))   
     if max(newline) gt 1 or min(newline) lt 0 then goto, tryscagain
     if(n_elements(newline) eq n_stars_per_clump[thisblock]) then $
     for s=0,n_stars_per_clump[thisblock]-1 do flagblock[3,bb[s]]=newline[s] $
     else goto,tryscagain
   endif
   print,''
; modify the FILENUMBERS if requested  ...
   print,'Do you wish to add/remove/change individual files y/N? > '
   d = get_kbrd(1) & wait,.5 
   if(d eq 'y' or d eq 'Y') then begin
      repeat begin
      print,'|Index Block Star N_Files SrcCal Dither  Type | Data Filenumbers...'
      for s=0,n_stars_per_clump[thisblock]-1 do begin
          w=where(fileblock[*,bb[s]] gt 0)
          print,flagblock[1,bb[s]],flagblock[2,bb[s]],s,flagblock[0,bb[s]],flagblock[3,bb[s]], $
                 flagblock[4,bb[s]],flagblock[5,bb[s]],fileblock[w,bb[s]],format="(3I5,5I7,100I5)"
      endfor
      print,'Input the star# of the line to be edited'
      read, instr, prompt='> '
      thistar=fix(instr)
      endrep until thistar ge 0 and thistar le n_stars_per_clump[thisblock]-1
      print,'Data Filenumbers are ...'
      w=where(fileblock[*,bb[thistar]] gt 0)
      print,fileblock[w,bb[thistar]],format="(90I5)"
      print,''
      print,'Paste modified numbers into command line'
      read, instr, prompt='> ' 
      newline=fix(string_vector(instr))   
      n_input=n_elements(newline)
      flagblock[0,bb[thistar]]=n_input   ; #number of files
      fileblock[*,bb[thistar]]=0
      fileblock[0:n_input-1,bb[thistar]]=newline
   endif
   
   ;add flags for calc_bispect, calibrate_v2_cp and binary_grid
   print,""
   print,'Do you wish to add or edit parameters for calc_bispect, calibrate_v2_cp and/or binary_grid y/N? > '
   d = get_kbrd(1) & wait,.5 
   if(d eq 'y' or d eq 'Y') then begin   
     print,'Input new parameters for calc_bispect if required, starting with a comma' 
     read, instr, prompt='>'    
     bispectflags[bb] = instr
     print,'Input new parameters for calibrate_v2_cp if required, starting with a comma' 
     read, instr, prompt='>'           
     calv2cpflags[bb] = instr
     print,'Input new parameters for binary_grid if required, starting with a comma' 
     read, instr, prompt='>'    
     bingridflags[bb] = instr
   endif
end


; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 4 %%%%%%%%%%%%%%%%%%%%%%%%%

;not working, disabled by changing case to '44'
;see note above and use delete then define new instead.

  '44': begin   ; overwrite a block 
   if n_clumps le 1 then begin
      print, '#### Number of defined blocks must be greater than 1.'
      print, '#### Please delete then add new.'
      goto, selectmode
   endif
   repeat begin
   print,'Input the Block# of the block to be overwritten. 0 -'+strcompress(string(n_clumps-1))
   read, instr, prompt='> ' 
   endrep until fix(instr) le n_clumps-1 and fix(instr) ge 0
   
   thisblock=fix(instr) 
   bb=where(flagblock[2,*] eq thisblock)

   ; Remove old defined block from arrays...
   not_bb=where(flagblock[2,*] ne thisblock)
     fileblock=fileblock[*,not_bb]
     flagblock=flagblock[*,not_bb]
     s_block=s_block[*,not_bb]
     n_block=n_block[*,not_bb]


   oldstarname=(starnames[bb])[0]
   olddatename=(datenames[bb])[0]
   starnames=starnames[not_bb]
   datenames=datenames[not_bb]

   starname=''
   Print,'Input star identifier for the new grouping, (default: ',oldstarname,')'
   read, starname, prompt=' > ' 
   if(starname eq '') then starname=oldstarname
   datename=''
   Print,'Input date identifier for the new grouping, (default: ',olddatename,')'
   read, datename, prompt=' > ' 
   if(datename eq '') then datename=olddatename
   repeat begin
   print,'Type the indexes of the stars in star grouping number ',strtrim(string(n_clumps),2)
   instr=''
   read, instr, prompt='> ' 
   thisclump=fix(string_vector(instr))   
   wc1 = [where(thisclump ge n_lines)]
   wc2 = [where(thisclump lt 0)]
   endrep until wc1 eq -1 and wc2 eq -1
   n_stars_thisclump=n_elements(thisclump)
   n_stars_per_clump[thisblock]=n_stars_thisclump
   for s=0, n_stars_thisclump-1 do begin 
     zvect=intarr(100)
     flagvect=intarr(10)
     nstarfiles=filenumlimits[1,thisclump[s]]-filenumlimits[0,thisclump[s]]+1
     thisfilevect=indgen(nstarfiles) + filenumlimits[0,thisclump[s]]
     zvect[0:nstarfiles-1]=thisfilevect  ; use array elements 10+ for the filenums
      ; Array elements 0-9 are reserved for information...
     flagvect[0]=nstarfiles    ; num files this star
     flagvect[1]=thisclump[s]  ; #index number through night
     flagvect[2]=thisblock     ; The Block number of this clump
     flagvect[3]=1             ; Default setting source/cal flag
     ; Try to figure out if it is a source or cal
     for t=0,nsrc_names-1 do $
         if( strpos(s_output[1,thisclump[s]],src_names[t]) ne -1) then flagvect[3]=0  ; Prolly a source
     if( strpos(s_output[1,thisclump[s]],'PSF') ne -1) then flagvect[3]=1             ; Prolly a CAL
     flagvect[4]=1             ; 0=nodither; 1=dither
     flagvect[5]=1             ; Observation Type:  1=normal src/cal;   2=faint companion hunt
     fileblock=[[fileblock],[zvect]] 
     flagblock=[[flagblock],[flagvect]]
     s_block  =[[s_block],[s_output[*,thisclump[s]]]]
     n_block  =[[n_block],[n_output[*,thisclump[s]]]]
     starnames = [starnames,starname]
     datenames = [datenames,datename]
   endfor
   end

; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 5 %%%%%%%%%%%%%%%%%%%%%%%%%
  '5': begin   ;  Delete block
   if n_clumps le 0 then begin
     print, "No Blocks Defined"
     goto, selectmode
   endif
   
   repeat begin
     print,'Do you really want to delete a block y/n? > '
     d = get_kbrd(1) & wait,.5 
   endrep until d eq 'y' or d eq 'Y' or d eq 'n' or d eq 'N'  
   if d eq 'n' or d eq 'N' then goto, selectmode

   ; if there is only 1 block defined, remove all defined variables
   if n_clumps eq 1 then begin
     ;delvar, fileblock, flagblock, s_block, n_block, datenames, starnames, skies, bispectflags, $
     ;        calv2cpflags, bingridflags, thisclump,n_stars_thisclump, n_stars_per_clump
     n_clumps = 0
     goto, selectmode
   endif
   
   repeat begin
   print,'Input the Block# of the block to be deleted. 0 -'+strcompress(string(n_clumps-1))
   read, instr, prompt='> ' 
   endrep until fix(instr) le n_clumps-1 and fix(instr) ge 0 
   thisblock=fix(instr) 
   notbb=where(flagblock[2,*] ne thisblock)

   fileblock    = fileblock[*,notbb]
   flagblock    = flagblock[*,notbb]
   s_block      = s_block[*,notbb]
   n_block      = n_block[*,notbb]
   datenames    = datenames[notbb]
   starnames    = starnames[notbb]
   skies        = skies[notbb]
   bispectflags = bispectflags[notbb]
   calv2cpflags = calv2cpflags[notbb]
   bingridflags = bingridflags[notbb]
   n_stars_per_clump=n_stars_per_clump[notbb]
   
   n_clumps = n_clumps-1
  end
   
   
; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 8 %%%%%%%%%%%%%%%%%%%%%%%%%
  '8': begin   ;  Save markup workings ...
   repeat begin
     print,'Save present markups y/n? > '
     d = get_kbrd(1) & wait,.5 
   endrep until d eq 'y' or d eq 'Y' or d eq 'n' or d eq 'N'
   
   if(d ne 'n' or d ne 'N') then begin
     save,n_clumps,skies,filenumlimits, fileblock,flagblock,s_block,n_block,starnames,datenames,n_stars_per_clump,$
          s_output,n_output,n_lines,new_src,filen,n_files,bispectflags,bispectflags,calv2cpflags,bingridflags,file,file=qb_name
     print,'Saved to file ',qb_name
   endif
   end

; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 9 %%%%%%%%%%%%%%%%%%%%%%%%%
  '9': begin   ; Done!
      start9:
      instr=''
      print,'Input which blocks you wish queued (a = all)' 
      read, instr, prompt='> '
      if(instr eq 'a' or instr eq 'A') then begin
        n_queuejobs=n_clumps
        queuejobs=indgen(n_clumps)
      endif else begin
        queuejobs=fix(string_vector(instr))
	if max(queuejobs) gt n_clumps-1 or min(queuejobs) lt 0 then goto, start9
        print, "processing the following blocks, " +string_vector(instr)
        n_queuejobs=n_elements(queuejobs)
      endelse 
      instr=''
      print,'Default perform analysis automatically, to output manual commands type m' 
      read, instr, prompt='> '
      if(instr eq 'm' or instr eq 'M') then manual=1 else manual=0
     goto,dothedata
   end

  'q': begin   ;  clean exit
    goto, last
  end

   else: goto,selectmode
endcase

goto,selectmode




; %%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis Starts Here %%%%%%%%%%%%%%%%%%%%%%%%%
dothedata:
if n_queuejobs le 0 then begin
    print, "No Jobs selected"
    goto, selectmode
endif

for job=0,n_queuejobs-1 do begin
  b=queuejobs[job]
  fbix=where(flagblock[2,*] eq b)  ; find the set belonging to block b
  bix=flagblock[1,fbix]            ; These are the actual indexes 
  
  ; Now we can find information about this particular configuration from the header data
  filtername=s_output(4,bix)
  maskname=s_output(6,bix)
  fn = where(filtername eq "empty")
  if fn[0] ne -1 then filtername[fn] = s_output(5,bix[fn])
  xydim=n_output(2,bix)
  xystr=strtrim(string(fix(xydim[0])),2)

  ; Check that the data bundled makes sense (same filter, same xydim)
  if( (size(uniq(filtername)))[0] ne 0  or (size(uniq(xydim)))[0] ne 0 ) then begin
    print,' ##############################################################################'
    print,' #####   WARNING - Incompatible filter or readout for Block ',strtrim(string(b),2),'             #####'
    print,' ##############################################################################'
    goto, selectmode
  endif
  ; Check that the data has a mask
  if( min(strpos(maskname,"Holes")) eq -1) then begin
    print,' ##############################################################################'
    print,' #####   WARNING - No mask in the beam with this data! Block ',strtrim(string(b),2),'            #####'
    print,' ##############################################################################'
    goto, selectmode
  endif
  ; Check for mask consistency
  if( (size(uniq(maskname)))[0] ne 0) then begin
    print,' ##############################################################################'
    print,' #####   WARNING - Incompatible mask for Block ',strtrim(string(b),2),'                          #####'
    print,' ##############################################################################'
    goto, selectmode
  endif

  theband=''
; add to here as required
  if( strpos(filtername[0],'J') ne -1 ) then theband='J' 
  if( strpos(filtername[0],'K') ne -1 ) then theband='K' 
  if( strpos(filtername[0],'H') ne -1 ) then theband='H' 
  if( strpos(filtername[0],'L') ne -1 ) then theband='L' 
  if( strpos(filtername[0],'IB_2.24') ne -1 ) then theband='IB224' 
  if( strpos(filtername[0],'NB_1.04') ne -1 ) then theband='NB104'
  if( strpos(filtername[0],'NB_1.28') ne -1 ) then theband='NB128'
  if( strpos(filtername[0],'NB_1.75') ne -1 ) then theband='NB175' 
  if( strpos(filtername[0],'NB_3.74') ne -1 ) then theband='NB374' 
  if( strpos(filtername[0],'NB_4.05') ne -1 ) then theband='NB405' 
 
  if theband eq '' then begin
    ; if this error occurs, the filter in the observation is not in the list above
    ; add it to the list as per the others and it will be resolved
    print, '##############################################################################'
    print, '#####   WARNING - Filter unrecognised by qball                           #####'
    print, '##############################################################################'
    goto, selectmode
  endif

  flatpath=flats_dir+'flat'+xystr+'_'+theband+'.idlvar'
  if file_test(flatpath) ne 1 then begin
    print, '##############################################################################'
    print, '#####   WARNING - Flat file '+flatpath+' not found'
    print, '##############################################################################'
    goto, selectmode
  endif
  ; Glue all the filenames into a vector and make up tsize vector
  frames=0 & tsize=0
  n_targets=n_elements(fbix)
  for f=0,n_targets-1 do begin 
     this_nfiles=flagblock(0,fbix[f])
     frames =  [frames,fileblock(0:this_nfiles-1,fbix[f]) ]
     tsize = [tsize,(fltarr(this_nfiles)+.1) * (f+1.)/10 * sign(flagblock[3,fbix[f]]-.5)]
  endfor
  frames =  frames[1:*]
  tsize =  tsize[1:*]

  ; Do we need to make up  a fake sky, or is it dithered?
  ditherstatus=flagblock[4,fbix] 
  if ditherstatus[0] eq 0 then begin
    skyfiles = *skies[fbix[0]].files
    skypath  = *skies[fbix[0]].path
  endif
  ;  COADD  = n_output(2,bix) 
  ;  MODE = n_output(3,bix) 
  ;  TINT = n_output(4,bix)
  ; OK I give up ... this finding an appropriate dark is just too much trouble.
  ; For the record, we need to pick one with coadds, cds, tint and xydim all right.
  ; This is way too annoying, and I am going to cop out like everyone else - PGT
    
  jobname=(starnames[fbix])[0]+'_'+(datenames[fbix])[0]  ; sinle name to identify this analysis

  ;expand prefix and updown, this is quick and rough, needs to be changed to allow for different prefixes
  prefix = replicate(prefix[0],n_elements(frames))
  updown = replicate(updown[0],n_elements(frames))
    

  if(manual eq 1) then begin     ; DO NOT do analysis ... just print out commands
      fstringi = '(a,"[",'+strcompress(string(n_elements(frames)),/rem)+'(I4,","),"]")'
      fstringf = '(a,"[",'+strcompress(string(n_elements(tsize)),/rem)+'(f5.2,","),"]")'
      print,'data_dir=',data_dir
      print,'jobname=',jobname
      print,'flatpath=',flatpath
      print,'frames=',frames,format=fstringi
      print,'tsize=',tsize,format=fstringf
      print,'extn=',extn
      if(ditherstatus[0] eq 0) then begin               ; NOT dithered data - skies (actually crap darks)
        print,'skypath =',skypath
        print,'skyfiles=',skyfiles
        print,'-----------'
        print,'qbe_conica, data_dir, jobname, flatpath, frames, skyfiles, tsize, prefix=prefix, extn=extn, $'
	print,'            updown=updown, horiz_destripe=horiz_destripe, discard_last = discard_last, ddir_sky=skypath'
      endif else begin 
        print,'-----------'
        print,'qbe_conica, data_dir, jobname, flatpath, frames, -1, tsize, prefix=prefix, $'
	print,'            extn=extn, updown=updown, horiz_destripe=horiz_destripe, discard_last = discard_last'
      endelse
      print,'calc_bispect, cubeinfo',jobname,'.idlvar , root_dir'+bispectflags[fbix[0]]
      print,'calibrate_v2_cp, cubeinfo',jobname,'.idlvar , root_dir'+calv2cpflags[fbix[0]]
      print,'restore, cubeinfo',jobname,'.idlvar'
      print,'-----------'
      ;goto, selectmode
    endif else begin

    ; OK here we actually run the main analysis 
    
      if file_test(original_dir+'/'+(starnames[fbix])[0],/directory) ne 1 then spawn, 'mkdir '+original_dir+'/'+(starnames[fbix])[0]
      cd,original_dir+'/'+(starnames[fbix])[0]
      if(ditherstatus[0] eq 0) then begin               ; NOT dithered data - skies (actually crap darks)
         qbe_conica, data_dir, jobname, flatpath, frames, skyfiles, tsize, prefix=prefix, extn=extn, $
	    	     updown=updown, horiz_destripe=horiz_destripe, discard_last = discard_last, ddir_sky=skypath
      endif else begin 
         qbe_conica, data_dir, jobname, flatpath, frames, -1, tsize, prefix=prefix, extn=extn, $
	    	     updown=updown, horiz_destripe=horiz_destripe, discard_last = discard_last
      endelse
      
      bispect = execute("calc_bispect, 'cubeinfo'+jobname+'.idlvar', root_dir"+bispectflags[fbix[0]])
      calv2cp = execute("calibrate_v2_cp, 'cubeinfo'+jobname+'.idlvar', root_dir"+calv2cpflags[fbix[0]])
      
      restore, 'cubeinfo'+jobname+'.idlvar'
      if((flagblock[5,fbix])[0] eq 2) then bingrid = execute('binary_grid, clog.primary_oifits'+bingridflags[fbix[0]]) else bingrid = -1

      cd,original_dir
    endelse

endfor
last:
cd,original_dir
end
