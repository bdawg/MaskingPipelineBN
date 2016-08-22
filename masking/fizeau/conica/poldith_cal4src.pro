;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;      poldith_cal4src
;
; PURPOSE:
;      Produces a cal4src matrix such that sources are only calibrated against         
;      calibrators in the same chip quadrant. Currently, this procedure only 
;      works for vertically polarized data that has been taken with successive 
;      left-right dither positions and it is assumed that there are equal numbers        
;      of source and calibrator data cubes. 
;
; CALLING SEQUENCE:
;      poldith_cal4src,nsrc,cal4src
;
; INPUTS:
;      nsrc - number of sources
;     
; OUTPUT:
;      cal4src - 
;
; HISTORY:
;      version 1.0  TME 11Jun09  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro poldith_cal4src,nsrc,cal4src


cal4src=fltarr(4*nsrc,4*nsrc)

ngroups=nsrc/2

for i=0,ngroups-1 do begin
  for j=0,3 do begin
    for k=j,nsrc*2-1,4 do begin
      cal4src[2*nsrc+k,4*i+j]=1
    endfor
  endfor
endfor

 

end
