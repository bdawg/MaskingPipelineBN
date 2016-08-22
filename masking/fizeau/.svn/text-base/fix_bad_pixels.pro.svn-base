;######## RELEASE CLIP HERE ######## 
function fix_bad_pixels,image,bad_pixels=bad_pixels


; bad_pixels is a list of bad pixels.

; ** NEW FEATURE **
;  if image pixels are NaN, then these will also be replace, and will
; display a warning message!
; JDM 14Mar99

if (keyword_set(image) eq 0) then begin
 print,'function fix_bad_pixels,im,bad_pixels=bad_pixels
 print,'    bad_pixels is a list of bad pixels.
return,-1
endif

; First thing to do is make a map of bad pixels for index file.
bad_mask=image
bad_mask(*)=0.0
bad_mask(bad_pixels)=1.0
nan_in=where( finite(image) eq 0,ct)
if (ct gt 0) then begin
   print,'** WARNING **'
   print,ct,' Pixels were NaN!! '
   print,' Correcting them via badpixel method..'
   bad_mask(nan_in)=1.0
   image(nan_in)=0.0
endif


; This will fill in bad pixels with average nearest neighbors
; If there are pixels which do not have any nearest neighbors, then it will
; increase the size of the neighborhood until all pixels do.  This is
; almost, but not QUITE the same at iterating with a constant size neighbor-
; hood.

; Initialize
new_bad_mask=bad_mask
newimage=image
num_still_bad=n_elements(bad_pixels)
still_bad=bad_pixels
neighborhood=1.0

while (num_still_bad gt 0) do begin
neighborhood=neighborhood+2 

norm=filter_image((1.0-bad_mask),smooth=neighborhood,/all)
still_bad=where(norm lt 1e-5,num_still_bad)
fix_these=where(new_bad_mask eq 1 and norm gt 0.0)
new_bad_mask(fix_these)=0.0
newimage(fix_these)=(filter_image(image*(1.-bad_mask),smooth=neighborhood,$
      /all))(fix_these)/norm(fix_these)
;stop

endwhile

;print,'To Fix all Bad Pixels, we used neighborhood diameter of ',$
;   neighborhood,' Pixels.'

return,newimage
end



