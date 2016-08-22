; 2001Sep19 JDM
; pro do_aperture_photometry, im, photometry
;
; This routine will take an image and do aperture photometry 
; based on 10 different hanning windows.
; This will assume sky subtraction has already been done.
;


pro photometry_nirc, im, photometry

; im is assumed to be centered at dimx/2, dimy/2 as it is
; in get_full and get_faint pow

aperture_photometry =[.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
num_ap=n_elements(aperture_photometry)
dimx=(size(im))(1)
dimy=(size(im))(2)

photometry=fltarr(num_ap)

;hans=fltarr(dimx,dimy,num_ap)
for i=0,num_ap-1 do begin
 diam = fix(sqrt(dimx*dimy)*aperture_photometry(i) )
 han0=hanning(diam,diam)
 han=fltarr(dimx,dimy)
 han(0:(diam-1),0:(diam-1))=han0
 locate_peak,han,x0,y0
 han=shift(han,-x0+dimx/2,-y0+dimy/2)
; hans(*,*,i) = han
 photometry(i)= total(im*han)
endfor
end

