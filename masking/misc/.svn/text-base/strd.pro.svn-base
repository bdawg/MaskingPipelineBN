PRO STRD, IM, HD, FILENAME        ;Read in SDAS image array and header array
;+
; NAME:
;	STRD
; PURPOSE:
;	Open an STSDAS file and read into an image array and header.  Combines
;	the functions of SXREAD and SXOPEN.  Can only be used on files 
;	without group parameters
;
; CALLING SEQUENCE:
;	STRD, im, hdr, filename  
;
; OPTIONAL INPUT:
;	FILENAME -  Character string giving the name of the SDAS file
;		to be read.  If omitted, then program will prompt 
;		for the file name.  If an extension is given, then
;		it must terminate in a 'h'.
;		A default extension of '.hhh' is assumed, if one is
;		not supplied.  VMS Version numbers are ignored, and the
;		most recent version is always used.
;
; OUTPUTS:
;	IM - array containing image data
;	HDR - string array containing header
;
; COMMON BLOCKS:
;	STCOMMN - Created by SXOPEN.  STRD uses STCOMMN to check
;		for an open unit, and to get image dimensions.          
;
; SIDE EFFECTS:
;	STSDAS image array and header are read into IM and HD
;	IF FILENAME is not supplied, then the program will check that
;	the image and header variable do not already contain data.
;
; RESTRICTIONS:
;	For use only on data without Groups!!  
;
; SYSTEM VARIABLES:
;	If !QUIET = 1 then program will not print the size of the image.
;
; PROCEDURE:
;	Program checks that STSDAS files exists and that IDL variables do
;	not already contain data, before calling SXOPEN and SXREAD to
;	read in SDAS data.
;
; MODIFICATION HISTORY:
;	Written W. Landsman, STI Corporation August 1986
;	Optional parameter "FILENAME" added November 1986
;-
 common stcommn,result,fname
 npar = N_params()
 err = string(7b) + 'STRD: ERROR - ' 
 warn = string(7b) + 'STRD: WARNING - '
 
 if npar eq 0 then begin
   print,'Syntax - strd,im,hd,[filename]
   return
 endif

if ( npar EQ 1 ) then begin    
   ans = ''
   print, warn, 'A name for the header array was not supplied.'
   read, 'Continue the procedure to read the image array only [YES]?', ANS
   if ( strmid(strupcase(ans),0,1) EQ 'N' ) then return    
endif     

if ( npar LT 3 ) and (N_elements(im) NE 0 ) then begin ;Already contain data?
   ans = ''   
   print, warn, 'Image array contains data that will be erased'
   read, 'Continue the procedure to read the image array [YES]? ',ANS
   if ( strmid( strupcase(ans),0,1 ) eq 'N') then return   
endif 

if npar EQ 3 then if ( N_elements( filename ) EQ 0) then $
       	print, err,'Third parameter must be character string' $
	else begin
		file_name = filename
	      	goto, FINDER
	endelse 

NAMER:

 file_name = ''          ;Get file name if not supplied
 read,'Enter name of SDAS data file (no quotes): ',FILE_NAME
 if file_name EQ '' then return    

FINDER: 
 fdecomp, file_name, disk, dir, name, ext, ver   

 if ver NE '' then $
    file_name = disk + dir + name + '.' + ext          ;No Versions allowed

 if ext EQ '' then file_name = file_name + '.hhh' $     ;Use default extension?
 else if strupcase( strmid(ext,2,1) ) NE 'H' then begin     
      	print,err, "SDAS file-name extensions must end with 'h'"
        goto, NAMER 
 endif

 find = findfile( file_name, COUNT = I )    ;Does file exist?
 if I LT 1 then begin
    print, err, 'Unable to find ' + SPEC_DIR( file_name )    
    if npar EQ 3 then return   
    print, 'Please re-enter name of file, or [RETURN] to exit'
    GOTO, namer  
 endif 

 for i = 1, 9 do begin               ;Find an open unit between 1 and 9
   test = fstat(i)
   if test.open eq 0 then begin
      unit = i
      goto, OPENER      
   endif
 endfor

OPENER: 
  sxopen, unit, file_name, hd   
 xsiz = result(10,unit)  & ysiz = result(11,unit) 

 if not !QUIET THEN $
   print, 'Now reading '+ strtrim(xsiz,2)+' by ' +strtrim(ysiz,2)+' array'

 im = sxread(unit) 

 close, unit    
 return
 end
