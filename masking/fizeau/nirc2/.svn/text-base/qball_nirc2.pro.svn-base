;;TO DO:
;;1) Space-separated file numbers in headtxt_output file.
;;2) Problems with automatic binary grid call

;+
; NAME:
;   QBALL_NIRC2
;
; PURPOSE:
;   Guide user through data markup for a whole data directory before 
;   passing the cubes through the analysis pipeline.
;
; DESCRIPTION:
;   QBALL_NIRC2 pulls useful information from all headers of a set of
;   data files and tries to produce a set of cubes clumping files
;   together when it finds something changes. 
;
;   The user has the option of specifying variables necessary for
;   the analysis (see below) with the qball_nirc2 call at the
;   command line. If these are not provided, the user will be
;   prompted to enter them manually once the qball_nirc2 call has
;   been made.
;
;   The user will then be guided through the data markup process. 
;   Once the markup for the data directory has been finished the
;   queued data is passed into the pipeline. The user also has 
;   the option of saving the markup as an idlvar file for later use.
;
; CALLING SEQUENCE:
;   qball_nirc2,[root_dir=root_dir, data_dir=data_dir, flats_dir=flats_dir, $
;                skyfiles=skyfiles, src_names=src_names, extn=extn,         $
;                headtxt_output=headtxt_output]
;
; INPUT:
;   (all inputs are optional - see 'optional keyword inputs')
;
; OUTPUT:
;
; OPTIONAL KEYWORD INPUTS:
;   root_dir - 
;   data_dir - directory containing the data to be analysed
;   flats_dir - directory containing the necessary flat files
;   skyfiles - integer array containing the frame numbers of skies
;   src_names - string array containing source names
;   extn - data files extension eg. '.fits.gz' 
;   headtxt_output - file name for the header that is produced
;   dont_ask - will accept any input you give it without asking you
;               if it's ok. Useful for preserving sanity
;-


pro qball_nirc2, root_dir=root_dir, data_dir=data_dir, flats_dir=flats_dir, skyfiles=skyfiles, $
                src_names=src_names, extn=extn, headtxt_output=headtxt_output,dont_ask=dont_ask


if (keyword_set(root_dir) eq 0 and keyword_set(data_dir) eq 0 and keyword_set(flats_dir) eq 0 $
    and keyword_set(skyfiles) eq 0 and keyword_set(src_names) eq 0 and keyword_set(extn) eq 0 $
    and keyword_set(headtxt_output) eq 0) then begin
  print,'No arguments passed to qball'
  print,'Would you like to use an existing markup (y) or create new (n)?' 
  d=get_kbrd(1) & wait,.5
  if (d eq 'y' or d eq 'Y') then goto,ready_to_go
endif

;; Now set the default root_dir if we can
defsysv, '!ROOT_DIR', exists=exists
if exists then root_dir=!ROOT_DIR
 if (keyword_set(root_dir) and keyword_set(data_dir) and keyword_set(flats_dir)) then begin
    if keyword_set(dont_ask) then d='y' else begin
       print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       print,'The following directory paths were passed to qball:'
       print,'root_dir = ',root_dir
       print,'data_dir = ',data_dir
       print,'flats_dir = ',flats_dir
       print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       print,'Accept these y/n?' 
       d = get_kbrd(1)          ; & wait,.5
    endelse
   if (d eq 'n' or d eq 'N') then begin
     goto,choose_dirs_manually
   endif else begin
     goto,skyfiles
   endelse
 endif else begin
   print,'One or more directory paths (root_dir,data_dir,flats_dir) have not been specified in the qball call'
   goto,choose_dirs_manually
 endelse
 choose_dirs_manually:
 print,'Would you like to type in the directory paths (t) or browse for them using a GUI (b)?'
 d = get_kbrd(1) & wait,.5
 tb_type_again:
 if (d eq 't' or d eq 'T' or d eq 'b' or d eq 'B') then begin
   if (d eq 't' or d eq 'T') then begin
     root_dir = ''
     data_dir = ''
     flats_dir = ''
     print,'Enter the full path for root_dir:'
     read,root_dir,prompt=' > '
     print,'Enter the full path for data_dir:'
     read,data_dir,prompt=' > '
     print,'Enter the full path for flats_dir:'
     read,flats_dir,prompt=' > '
     goto,skyfiles
   endif
   if (d eq 'b' or d eq 'B') then begin
     print,''
     print,'Select the root_dir'
     root_dir = dialog_pickfile()
     print,'  you selected the root_dir: ',root_dir
     print,''
     print,'Select the data_dir'
     data_dir = dialog_pickfile()
     print,'  you selected the data_dir: ',data_dir
     print,'Select the flats_dir'
     flats_dir = dialog_pickfile()
     print,'  you selected the flats_dir: ',flats_dir
     goto,skyfiles
   endif
 endif else begin
   print,'Please type either t to type or b to browse'
   d = get_kbrd(1) & wait,.5
   goto,tb_type_again
 endelse


 skyfiles:
 if keyword_set(skyfiles) then begin
    if keyword_set(dont_ask) then d='y' else begin
       print,''
       print,'You have entered the following frames as sky files in the qball call:'
       print,skyfiles
       print,'Would you like to accept these y/n?'
       d = get_kbrd(1)
    endelse
    if (d eq 'y' or d eq 'Y') then begin
       goto,src_names
    endif 
 endif
 skyfiles = ''
 display_skyfiles = strarr(2)
 add_more_skyfiles: 
 these_skyfiles = ''
 read,these_skyfiles,prompt='Enter the first and last sky frames for the block of skies (separated by spaces) > '
 these_skyfiles = fix(string_vector(these_skyfiles))
 display_skyfiles = [[display_skyfiles],[these_skyfiles]]
 skyfiles = [skyfiles,these_skyfiles[0]+indgen(these_skyfiles[1]-these_skyfiles[0]+1)]
 if keyword_set(dont_ask) then d='n' else begin
    print,'Would you like to add more skyfiles y/n?'
    d = get_kbrd(1) & wait,.5
 endelse
 if (d eq 'y' or d eq 'Y') then begin
    goto,add_more_skyfiles
 endif else begin
   print,'You have entered the following skyfiles: '
   print,'   First     Last'
   print,display_skyfiles[*,1:*]
   if keyword_set(dont_ask) then d='y' else begin
      print,'Would you like to accept these (y) or start over (n)?'
      d = get_kbrd(1) & wait,.5
   endelse
   if (d eq 'n' or d eq 'N') then begin
      skyfiles = skyfiles[1:*]  ; remove the meaningless zero from the first entry
      goto,skyfiles
   endif else begin
      goto,src_names
   endelse 
endelse


 src_names:
 if keyword_set(src_names) then begin
    if keyword_set(dont_ask) then d='y' else begin
       print,''
       print,'You have entered the following source names in the qball call:'
       print,src_names
       print,'Would you like to accept these y/n?'
       d = get_kbrd(1) & wait,.5
    endelse
    if (d eq 'y' or d eq 'Y') then begin
       goto,extn
    endif 
 endif
 src_names = ''
 add_more_src_names: 
 these_src_names = ''
 read,these_src_names,prompt='Enter source names (separated by spaces) > '
 these_src_names = string_vector(these_src_names)
 src_names = [src_names,these_src_names]
 if keyword_set(dont_ask) then d='n' else begin
    print,'Would you like to add more src_names y/n?'
    d = get_kbrd(1) & wait,.5
 endelse
 if (d eq 'y' or d eq 'Y') then begin
    goto,add_more_src_names
 endif else begin
   print,'You have entered the following src_names:'
   print,src_names
   print,'Would you like to accept these (y) or start over (n)?'
   d = get_kbrd(1) & wait,.5
   if (d eq 'n' or d eq 'N') then begin
     src_names = ''
     goto,src_names
   endif else begin
     src_names = src_names[1:*]  ; remove the meaningless first empty entry
     goto,extn
   endelse 
 endelse


 extn:
 nsrc_names=n_elements(src_names)
 if keyword_set(extn) then begin
    if keyword_set(dont_ask) then d='y' else begin
       print,''
       print,'You have entered the following extn for data files in the qball call:'
       print,extn
       print,'Would you like to accept this extn y/n?'
       d = get_kbrd(1) & wait,.5
    endelse
    if (d eq 'y' or d eq 'Y') then begin
       goto,headtxt_output
    endif 
 endif
 extn = ''
 read,extn,prompt='Enter the full extn to be used for data files > '


 headtxt_output:
 if keyword_set(headtxt_output) then begin
    if keyword_set(dont_ask) then d='y' else begin
       print,''
       print,'You have entered the following headtxt_output in the qball call:'
       print,headtxt_output
       print,'Would you like to accept this headtxt_output y/n?'
       d = get_kbrd(1) & wait,.5
    endelse
    if (d eq 'y' or d eq 'Y') then begin
       goto,ready_to_go
    endif 
 endif

 headtxt_output=''
 read,headtxt_output,prompt='Enter the file name to be used for headtxt_output > ' 

 ; %%%%%%%%%%%%%%%%%%%%
 ; AUTOMATIC BELOW HERE
 ; %%%%%%%%%%%%%%%%%%%%

 ready_to_go:

 print,''
 print,'Now guiding you through the data markup...'
 print,''

 ; Find the current working directory
 CD, CURRENT=original_dir

 ; Housekeeping to save/restore previous qball markups ...
 qb_name=''
 read, qb_name, prompt='Name of this data markup (e.g. qb_nirc2_090531.idlvar) > ' 
 print,'Restore previous y/N? > '
 d = get_kbrd(1) & wait,.5 
 if(d eq 'y' or d eq 'Y') then begin
    restore,qb_name
    goto,skipheadread
 endif

 ; Things to grab from header ...
 s_par=['FILTER','OBJECT','TARGNAME','SLITNAME','CAMNAME','PMSNAME','ROTMODE','UTC']
 n_par=['NAXIS1','NAXIS2','COADDS','SAMPMODE','ITIME','EL','PARANG','ROTPPOSN','RA','DEC']

 if (keyword_set(data_dir) eq 0) then data_dir="./"
 file=findfile(data_dir+'n*'+extn)
 if file[0] eq "" then begin
     print, "No files found"
     goto, last
 endif
 n_files=n_elements(file)
 spars=n_elements(s_par)
 npars=n_elements(n_par)
 s_output=strarr(spars,n_files)
 n_output=fltarr(npars,n_files)


 ; Read through all the headers, grab relevant info...
 print, strcompress(string(n_files)), " files found. Reading..."
 new_src=[0]
 for i=0,n_files-1 do begin
   print,i,file(i)
    filename=file(i)
    head=headfits(filename)
   for j=0,n_elements(s_par)-1 do begin
    s_output(j,i)=sxpar(head,s_par(j))
   endfor
   for j=0,n_elements(n_par)-1 do begin
    n_output(j,i)=sxpar(head,n_par(j))
   endfor
 ; work out if the present observation is a new target/color/nod
 ; or simply another data file same as before ...
   if(i gt 0) then                                $
     if( abs(n_output(8,i)-run_ra) gt 1. or       $   ; change in pointing?
         abs(n_output(9,i)-run_dec) gt .01 or     $   ;
         strjoin(s_output(0:6,i)) ne run_strings or $   ; change in config (strings)?
         strjoin(string(n_output(0:4,i))) ne run_ints ) $ ; change in config (int)?
       then new_src = [new_src,i]
   run_ra=n_output(8,i)
   run_dec=n_output(9,i)
   run_strings=strjoin(s_output(0:6,i))
   run_ints=strjoin(string(n_output(0:4,i)))
 endfor   

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
 print,' IX  File#           FILTER      OBJECT     TARGET      SLIT CAMERA   PUPIL     ROTMODE         UTC     AX1  AX2 COAD MODE TINT  EL   PARANG ROTPPOSN'
 for i=0,n_lines-1 do begin
   pos=strpos(filen(i),'.fits')
   file_num=fix(strmid(filen(i),pos-4,4))
   filenumlimits[0,i]=file_num
   if(where(new_src eq new_src(i)+1) eq -1 and i lt n_lines-1) then begin
     pos2=strpos(file(new_src(i+1)-1),'.fits')
     file_interval='-'+ strmid(file(new_src(i+1)-1),pos-2,2)
     filenumlimits[1,i]=fix(strmid(file(new_src(i+1)-1),pos-4,4))
   endif else file_interval=''  
   if(i eq n_lines-1 and new_src(i) ne n_files-1) then file_interval='-'+ strmid(file(n_files-1),pos-2,2)
   print,i,file_num,file_interval,s_output(*,i),n_output(0:7,i),format="(2I4,A3,' ',A17,2(A12),A7,A7,A9,A14,A12,3(I5),I2,F8.3,F5.1,2F8.2)"
 endfor
 if(filenumlimits[1,n_lines-1] eq 0) then filenumlimits[1,n_lines-1]=fix(strmid(file(n_files-1),pos-4,4))

 ; Now write to output file
 openw,2,headtxt_output
 printf,2,' IX  File#           FILTER      OBJECT     TARGET      SLIT CAMERA   PUPIL     ROTMODE         UTC     AX1  AX2 COAD MODE TINT  EL   PARANG ROTPPOSN'
 for i=0,n_lines-1 do begin
   pos=strpos(filen(i),'.fits')
   file_num=fix(strmid(filen(i),pos-4,4))
   if(where(new_src eq new_src(i)+1) eq -1 and i lt n_lines-1) then begin
     pos2=strpos(file(new_src(i+1)-1),'.fits')
     file_interval='-'+ strmid(file(new_src(i+1)-1),pos-2,2)
     filenumlimits[1,i]=fix(strmid(file(new_src(i+1)-1),pos-4,4))
   endif else file_interval=''  
   if(i eq n_lines-1 and new_src(i) ne n_files-1) then file_interval='-'+ strmid(file(n_files-1),pos-2,2)
   printf,2,i,file_num,file_interval,s_output(*,i),n_output(0:7,i),format="(2I4,A3,' ',A17,2(A12),A7,A7,A9,A14,A12,3(I5),I2,F8.3,F5.1,2F8.2)"
 endfor
 close,2



 n_clumps=0

 skipheadread:

 selectmode:
 print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 print,'Select option:'
 print,'0 = Print original header strip'
 print,'1 = Define a new data Block'
 print,'2 = Show presently defined Blocks'
 print,'3 = Edit a Block'
 print,'4 = Overwrite a Block'
 print,'8 = Save present state of data markups'
 print,'9 = Finished with data tagging - queue analysis'
 print,'q = quit'
 print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 k = get_kbrd(1) & wait,.5
 case k of


 ; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 0 %%%%%%%%%%%%%%%%%%%%%%%%%
   '0': begin   ; Print stuff out...
 print,' IX Block File#           FILTER      OBJECT     TARGET      SLIT CAMERA   PUPIL     ROTMODE         UTC     AX1  AX2 COAD MODE TINT  EL   PARANG ROTPPOSN'
 for i=0,n_lines-1 do begin
   pos=strpos(filen(i),'.fits')
   file_num=fix(strmid(filen(i),pos-4,4))
   if(where(new_src eq new_src(i)+1) eq -1 and i lt n_lines-1) then begin
     pos2=strpos(file(new_src(i+1)-1),'.fits')
     file_interval='-'+ strmid(file(new_src(i+1)-1),pos-2,2)
   endif else file_interval=''  
   blockstring=''
   bw=where(flagblock[1,*] eq i) 
     if(bw[0] ne -1) then begin
       if(n_elements(bw) eq 1) then blockstring=strtrim(string(flagblock[2,bw]),2) $
         else print,'* NOTE - this star set in mutiple blocks!'
       bw=bw[0]
     endif
   if(i eq n_lines-1 and new_src(i) ne n_files-1) then file_interval='-'+ strmid(file(n_files-1),pos-2,2)
   print,i,blockstring,file_num,file_interval,s_output(*,i),n_output(0:7,i),format="(I4,A4,I4,A3,' ',A17,2(A12),A7,A7,A9,A14,A12,3(I5),I2,F8.3,F5.1,2F8.2)"
 endfor
 end


 ; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 1 %%%%%%%%%%%%%%%%%%%%%%%%%
   '1': begin   ; define a new group of files to calibrate together

    starname=''
    datename=''
    read, starname, prompt='Star identifier/working dir for the new grouping (e.g. EtaSer) > ' 
    read, datename, prompt='Date identifier for the new grouping (e.g. 25Jun09) > ' 
    repeat begin
    print,'Type the indexes of the stars in star grouping number (separated by spaces) ',strtrim(string(n_clumps),2)
    instr=''
    read, instr, prompt='> ' 
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
      for t=0,nsrc_names-1 do begin $
          if( strpos(s_output[2,thisclump[s]],src_names[t]) ne -1) then flagvect[3]=0  ; Prolly a source
          print,s_output[2,thisclump[s]],src_names[t],flagvect[3]
       endfor
      if keyword_set(skyfiles) then flagvect[4]=0 else flagvect[4]=1 ; 0=nodither; 1=dither
      flagvect[5]=1             ; Observation Type:  1=normal src/cal;   2=faint companion hunt
      if(s eq 0 and n_clumps eq 0) then begin
           fileblock=zvect 
           flagblock=flagvect
           s_block  =s_output[*,thisclump[s]]
           n_block  =n_output[*,thisclump[s]]
           datenames = datename
           starnames = starname
      endif else begin
           fileblock=[[fileblock],[zvect]] 
           flagblock=[[flagblock],[flagvect]]
           s_block  =[[s_block],[s_output[*,thisclump[s]]]]
           n_block  =[[n_block],[n_output[*,thisclump[s]]]]
           datenames = [datenames,datename]
           starnames = [starnames,starname]
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
    print,'Input the Block# of the block to be edited'
    read, instr, prompt='> ' 
    thisblock=fix(instr) 
    endrep until thisblock ge 0 and thisblock lt n_clumps
    bb=where(flagblock[2,*] eq thisblock)
   
 ; modify the DITHER and TYPE flags if requested  ...
    print,'This block: Dither, Type = ',strtrim(string(flagblock[4:5 ,bb[0]]),2),' : modify y/N? > '
    d = get_kbrd(1) & wait,.5 
    if(d eq 'y' or d eq 'Y') then begin
      repeat begin
      print,'Input new values for Dither and Type e.g. > 1 2 ' 
      read, instr, prompt='>' 
      newline=fix(string_vector(instr))  
      endcond = 0
      if n_elements(newline) eq 1 and newline[0] eq 0 then endcond = 1 else $
      if n_elements(newline) eq 2 and newline[0] eq 1 and (newline[1] eq q or newline[1] eq 2) then endcond = 1
      endrep until endcond eq 1
        flagblock[4,bb]=newline[0]
        if newline[0] eq 1 then flagblock[5,bb]=newline[1] else flagblock[5,bb]=1
    endif
    print,''
 ; modify the SRC/CAL flag if requested  ...
    print,'This block Src/Cal sequence = ',strtrim(string(reform(flagblock[3,bb])),2)
    print,'Modify y/N? > '
    d = get_kbrd(1) & wait,.5 
    if(d eq 'y' or d eq 'Y') then begin
      tryscagain:
      print,'Input new values for Src/Cal flag 0=Src, 1=Cal e.g. >  0 1 0 1 0 1 ' 
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
 end


 ; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 4 %%%%%%%%%%%%%%%%%%%%%%%%%
   '4': begin   ; overwrite a block 
    if n_clumps le 1 then begin
       print, '#### Number of defined blocks must be greater than 1 ####'
       goto, selectmode
    endif
    repeat begin
    print,'Input the Block# of the block to be overwritten'
    read, instr, prompt='> ' 
    thisblock=fix(instr) 
    endrep until fix(instr) le n_clumps-1 and fix(instr) ge 0
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
    wc = [[where(thisclump ge n_lines)],[where(thisclump lt 0)]]
    endrep until wc eq -1
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
      flagvect[2]=thisblock      ; The Block number of this clump
      flagvect[3]=1             ; Default setting source/cal flag
      ; Try to figure out if it is a source or cal
      for t=0,nsrc_names-1 do $
          if( strpos(s_output[2,thisclump[s]],src_names[t]) ne -1) then flagvect[3]=0  ; Prolly a source
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


 ; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 8 %%%%%%%%%%%%%%%%%%%%%%%%%
   '8': begin   ;  Save markup workings ...
    print,'Save present markups Y/n? > '
    d = get_kbrd(1) & wait,.5 
    if(d ne 'n' or d ne 'N') then begin
      save,root_dir,data_dir,flats_dir,skyfiles,src_names,extn,headtxt_output, $
         n_clumps,fileblock,flagblock,s_block,n_block,starnames,datenames,n_stars_per_clump,$
         s_output,n_output,n_lines,new_src,filen,n_files,file,file=qb_name
      print,'Saved to file ',qb_name
    endif
    end

 ; %%%%%%%%%%%%%%%%%%%%%%%%% CASE 9 %%%%%%%%%%%%%%%%%%%%%%%%%
   '9': begin   ; Done!
       instr=''
       print,'Input which blocks you wish queued (a = all)' 
       read, instr, prompt='> '
       if(instr eq 'a' or instr eq 'A') then begin
         n_queuejobs=n_clumps
         queuejobs=indgen(n_clumps)
       endif else begin
         queuejobs=fix(string_vector(instr))
         n_queuejobs=n_elements(queuejobs)
       endelse 
       instr=''
       print,'Default perform analysis automatically, to output manual commands type m' 
       read, instr, prompt='> '
       if(instr eq 'm' or instr eq 'M') then manual=1 else manual=0
      goto,dothedata
    end

   'q':begin
       goto, last
    end

    else: goto,selectmode


 endcase


 goto,selectmode




 ; %%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis Starts Here %%%%%%%%%%%%%%%%%%%%%%%%%
 dothedata:



 for job=0,n_queuejobs-1 do begin
   b=queuejobs[job]
   fbix=where(flagblock[2,*] eq b)  ; find the set belonging to block b
   bix=flagblock[1,fbix]            ; These are the actual indexes 
  
   ; Now we can find information about this particular configuration from the header data
   filtername=s_output(0,bix)
   xydim=n_output(0,bix)
   xystr=strtrim(string(fix(xydim[0])),2)

   ; Check that the data bundled makes sense (same filter, same xydim)
   if( (size(uniq(filtername)))[0] ne 0  or (size(uniq(xydim)))[0] ne 0 ) then begin
     print,' ##############################################################################'
     print,' ##### WARNING - Incompatible filter or readout for Block',strtrim(string(b),2),' #####'
     print,' ##############################################################################'
   endif
   ; Check that the data has a mask 
   if( strpos(filtername[0],'18h') eq -1  and strpos(filtername[0],'9h') eq -1 ) then begin
     print,' ##############################################################################'
     print,' ##### WARNING - No mask in the beam this data!',strtrim(string(b),2),' #####'
     print,' ##############################################################################'
   endif

   thefilter=''
   if( strpos(filtername[0],'J') ne -1 ) then theband='J'
   if( strpos(filtername[0],'K') ne -1 ) then theband='K'
   if( strpos(filtername[0],'H') ne -1 ) then theband='H'
   if( strpos(filtername[0],'L') ne -1 ) then theband='L'

   flatpath=flats_dir+'flat_'+theband+'_'+xystr+'x'+xystr+'.idlvar'
 
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

   ;  COADD  = n_output(2,bix) 
   ;  MODE = n_output(3,bix) 
   ;  TINT = n_output(4,bix)
   ; OK I give up ... this finding an appropriate dark is just too much trouble.
   ; For the record, we need to pick one with coadds, cds, tint and xydim all right.
   ; This is way too annoying, and I am going to cop out like everyone else - PGT
    
   jobname=(starnames[fbix])[0]+'_'+(datenames[fbix])[0]  ; sinle name to identify this analysis

   if(manual eq 1) then begin     ; DO NOT do analysis ... just print out commands
       print,'data_dir="',data_dir,'"'
       print,'root_dir="',root_dir,'"'
       print,'jobname="',jobname,'"'
       print,'flatpath="',flatpath,'"'
       print,'frames=',frames,format='(a,"[",'+strcompress(string(n_elements(frames)),/remove_all)+'(I4,","),"]")'
       print,'tsize=',tsize,format='(a,"[",'+strcompress(string(n_elements(tsize)),/remove_all)+'(f5.2,","),"]")'
       print,'extn="',extn,'"'
       if(ditherstatus[0] eq 0) then begin               ; NOT dithered data - skies (actually crap darks)
         print,'skyfiles=',skyfiles,format='(a,"[",'+strcompress(string(n_elements(frames)),/remove_all)+'(I4,","),"]")'
         print,';-----------'
         print,'qbe_nirc2, data_dir, jobname, flatpath, skyfiles, fr=frames, ts=tsize, extn=extn'
       endif else begin 
         print,';-----------'
         print,'qbe_nirc2, data_dir, jobname, flatpath, -1, fr=frames, ts=tsize, extn=extn'
       endelse
       print,'calc_bispect, "cubeinfo',jobname,'.idlvar" , root_dir=root_dir'
       print,'calibrate_v2_cp, "cubeinfo',jobname,'.idlvar" , root_dir=root_dir'
       print,'restore, "cubeinfo',jobname,'.idlvar"'
       print,';-----------'
       print,'end'
     endif else begin

     ; OK here we actually run the main analysis 
       spawn,'cd '+original_dir+'/'+(starnames[fbix])[0],dir_there
       if (keyword_set(dir_there) eq 0) then $
         spawn,'mkdir '+original_dir+'/'+(starnames[fbix])[0]
       analysis_dir = original_dir+'/'+(starnames[fbix])[0]
       print,'Would you like the analysis output saved in the current directory (c)'
       print,'or would you like it saved in '+analysis_dir+' (s)'
       print,'or would you like to specify another directory (n)?'
       d = get_kbrd(1) & wait,.5
       if (d eq 'c' or d eq 'C') then begin
         spawn,'pwd',analysis_dir
       endif
       if (d eq 'n' or d eq 'N') then begin
         print,'Enter the directory path for the analysis output to be sent to:'
         read,analysis_dir,prompt=' > '
         spawn,'cd '+analysis_dir,dir_there
         if (keyword_set(dir_there) eq 0) then $
           spawn,'mkdir '+analysis_dir
       endif
       cd,analysis_dir
       if(ditherstatus[0] eq 0) then begin               ; NOT dithered data - skies (actually crap darks)
          qbe_nirc2, data_dir, jobname, flatpath, skyfiles, fr=frames, ts=tsize, extn=extn
       endif else begin 
          qbe_nirc2, data_dir, jobname, flatpath, -1,       fr=frames, ts=tsize, extn=extn
       endelse
 
       calc_bispect, 'cubeinfo'+jobname+'.idlvar', root_dir
       calibrate_v2_cp, 'cubeinfo'+jobname+'.idlvar', root_dir
       restore, 'cubeinfo'+jobname+'.idlvar'

       if(flagblock[5,fbix] eq 2) then $
        binary_grid, clog.primary_oifits

       cd,original_dir
     endelse

 endfor
 last:
 cd,original_dir
end
