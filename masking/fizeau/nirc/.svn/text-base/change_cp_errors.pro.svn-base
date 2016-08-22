function change_cp_errors,old_cp_errors,bad_flag=bad_flag,bad_num=bad_num

;  JDM early 1997
;
;  new_cp_errors=change_cp_errors(old_cp_errors,bad_num=bad_num,
;   bad_flag=50)  (in degrees!)
;     	
;
; 
;    The IDL routines keep track of the "error" in a phase measurement in
;    a very specific mathematical sense.  It is tan of the ratio between
;    the stdev of a quantity and the quantity itself.  Hence perfectly
;    meaningless quantity will have an "error" of 90 degrees,
;    the angle whose tangent is infinity.  
;    
;
;     When these angles are put into the VLBI package, we should
;     increase the ERRORS to more consistent for fitting routines
;     which are looking at chi-squared statistics.  
;     This actually does not require very big changes, but changes
;     which should be made.  In addition, the BAD_FLAG optional input
;     variable will tell this routine to invalidate all closure phases
;     with IDL "ERROR" higher than this amount (in degrees).  
;     This might be beneficial, since even this routine say the
;     error on a given quantity is 120 degrees, even for SNR=0.
;	Also the variable bad_num is returned with the number of baselines which
; 	were invalidated due to the bad_flag
;
; PGT 05Nov98 - small change to bad_flag operation.
; JDM 31Mar01 - made this routine independent of the change_cp_errors.idlvar file by
; 	        embedding the monte carlo results in the code.  This was done to 
;		ease portability of code to different machines.
;
if (keyword_set(old_cp_errors) eq 0) then begin
print,'new_cp_errors=change_cp_errors(old_cp_errors,bad_num=bad_num,$
print,'      bad_flag=45)  (in degrees!)
return,0
endif
cp_errors=old_cp_errors
if (keyword_set(bad_flag) eq 0) then bad_flag=180.

goto,newmethod      ;  Most recent ( and hopefully accurate)
		    ;  method for performing this routine's task is
		    ;  found at the end. JDM

angles=findgen(90)+.5
weights=tan(angles/!radeg)

;angles_array=fltarr(90,360)
new_array=fltarr(90)

for i=0,89 do begin
   rp=1.0+weights(i)*cos(findgen(360)/!radeg)
   ip=weights(i)*sin(findgen(360)/!radeg)
  
   ri2at,rp,ip,a,t
   t=mod360(t)
   new_array(i)=stdev(t)
endfor
;_____________________NEW METHOD HERE___________________
; Based on Monte Carlo Study.   April 13, 1997. JDM
newmethod:

;new_array=fltarr(180)
;angles=fltarr(180)
;restore,root_dir+"Speckle/change_cp_errors.idlvar'
; The raw data is in Speckle/Closure_Phase_Work/change_cp_errors.idlvar, but I will
; use lower-resolution, smoothed version:
angles=[0.00, 3.07, 6.17, 9.31, 12.54, 16.07, 20.24, 25.19, 30.53, 35.52, 39.10, 41.6810, $
      43.08, 43.94, 44.46, 44.65, 44.65, 44.81, 45.54]
new_array =[0.00, 3.19, 6.32, 9.56, 12.87, 16.51, 21.27, 28.51, 39.10, 51.34, 63.30,$
      74.30, 83.09, 90.08, 94.56, 98.09, 100.50, 102.432, 105.460]



index=where(cp_errors ge bad_flag ,bad_num)
if (bad_num gt 0) then cp_errors(index)=-1.
index=where(cp_errors ne -1,ct)
if (ct gt 0) then $
cp_errors(index)=interpol(new_array,angles,cp_errors(index))
return,cp_errors
end




