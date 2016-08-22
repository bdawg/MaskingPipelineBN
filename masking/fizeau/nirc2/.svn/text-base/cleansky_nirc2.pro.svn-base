; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    cleansky_nirc2.pro                          %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program is intended for use with NIRC2 speckle data.
;     It does steps beginning with raw data and outputting a cleaned cube
;
; Input Variables:
;     skynum  : filennum of the sky frames.
;  bad_pixels : 
;    datadir  : data directory 
;    root_dir : directory to look for templates (duvet fix)
; Output:
;    supersky : cleaned, sky data frame
;  headinfo   : structure of info stripped from header
; Options:
;    /destripe:     Dstripes the image
;		    
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; History : Written                                         PGT  26Feb04
;           Extensive changes/cleaning up. 
;           Only works on a single data cube at a time


pro cleansky_nirc2,skynum,bad_pixels=bad_pixels,	$
        datadir=datadir,root_dir=root_dir,supersky=supersky,  extn = extn

cleanplot,/silent    ; Intialize Display Driver
if (keyword_set(datadir) eq 0) then datadir="./"
;if (keyword_set(gain) eq 0) then gain=1.0
if (keyword_set(bad_pixels) eq 0) then bad_pixels=0
if (keyword_set(root_dir) eq 0) then root_dir='~/code/'
if (keyword_set(extn) eq 0) then extn='.fits'

;______
; Initialization
;______

; %%%%%  Read in sky frames & Create Supersky
filename=datadir+"n"+string(skynum,format="(I4.4)")+extn
skies=float(reform(readfits(filename,skyhead)))
info=size(skies)
if(info[0] eq 2) then begin
   nsky=1 
   supersky=skies
endif else begin
   nsky=info(3)
   supersky=total(skies,3)/nsky
endelse

; %%%%%  Remove Bad Pixels 
supersky=fix_bad_pixels(supersky,bad=bad_pixels)
; %%%%% No Gain correction here:
; supersky=supersky/gain
print,'Supersky made up of: ',nsky,' Sky Frames'

end
