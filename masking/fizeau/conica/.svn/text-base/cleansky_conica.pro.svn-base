;######## RELEASE CLIP HERE ######## 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    cleansky_conica.pro                          %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program is intended for use with conica speckle data.
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


pro cleansky_conica,skynum,bad_pixels=bad_pixels,prefix=prefix,extn=extn,	$
        datadir=datadir,root_dir=root_dir,supersky=supersky

cleanplot,/silent    ; Intialize Display Driver
if (keyword_set(datadir) eq 0) then datadir="./"
;if (keyword_set(gain) eq 0) then gain=1.0
if (keyword_set(bad_pixels) eq 0) then bad_pixels=0
if (keyword_set(root_dir) eq 0) then root_dir='~/code/'
if (keyword_set(prefix) eq 0) then prefix =  'NACO_IMG_SCI'
if (keyword_set(extn) eq 0) then extn =  '_DIT.fits'

;______
; Initialization
;______

; %%%%%  Read in sky frames & Create Supersky
filename=datadir+prefix+string(skynum,format="(I4.4)")+extn
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
