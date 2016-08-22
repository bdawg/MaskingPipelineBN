; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; This program generates useful information for masks
; Takes as input the (x,y) locations of the holes in mm on primary
; Currently, gives (x,y) pupil postions in meters
;                  (x,y) locations of Fourier Plane Coverage (in meters)
; Lengths of baselines in a big vector (in meters)
; Orientations of baselines in a big vector (in radians)

pro mask_info,xy_in,xy_coords,F_cov,b_lengths,b_angles

s=size(xy_in)
n_holes=s(1)

xy_coords=xy_in /1000  

;Now make up array F_cov of expected fourier coverage
F_cov=fltarr(2,n_holes,n_holes)
for i=0, n_holes-1 do begin
  for j=i+1, n_holes-1 do begin
    F_cov(0,i,j)=( xy_coords(i,0)-xy_coords(j,0) ) 
    F_cov(0,j,i)=-F_cov(0,i,j)
    F_cov(1,i,j)=( xy_coords(i,1)-xy_coords(j,1) )
    F_cov(1,j,i)=-F_cov(1,i,j)
  endfor
endfor
 
;now work out baselines
b_lengths = reform(sqrt( F_cov(0,*,*)^2 + F_cov(1,*,*)^2 ))

; Now do Angles
  b_angles  = reform(atan( F_cov(1,*,*) / (F_cov(0,*,*)+1e-10) ))
  b_angles(where( F_cov(0,*,*) lt 0))=  $
       b_angles(where(F_cov(0,*,*) lt 0)) + 3.1415926

end

