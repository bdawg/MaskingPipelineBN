; keckmrg_input.pro
;    
;  JDM  Written.. February 1997
;  JDM  Added /append   April 1997
;  JDM  Put in zero phases for short baselines first.	Oct 1998
;  PGT  Added ERROR_HISTOGRAM and clp error cutoff stuff Nov 1998
;  PGT  modified ERROR_HISTOGRAM to opereate before calc_phases
;  JDM  modified to find minimal set of baselines based on 
;	'bad' closure triangles (keyword: purge_baselines) Nov 29, 1998
;  PGT  Put /swap_if_little_endian on open statments .... (everything SUN format)
;  JDM  renamed keckmrg_input (formerly johnmrg_input) Mar 31, 2001
;
; This procedure is designed to create an ASCII file for JOHNMRG to read.
; JOHNMRG (a fortran routine) will create the BINARY MRG file for use 
; with the CALTECH VLBI package.
;
; This routine will output the following format:
;     NSTATS    (INT)   Number of Stations; limit 99 stations
;     WAV       (REAL)  Wavelength (meters)
;     NBAS      (INT)   Number of Baselines
;     NTR       (INT)   Number of Closure Trianges
;     (X(I),I=1,NSTATS) (REAL)  X-coord of Station (meters on the primary)
;     (Y(I),I=1,NSTATS) (REAL)  Y-coord of Station (meters on the primary)
;     (BAS[I], I=1,NBAS)(INT)   Caltech Notation
;                               100*(station1)+(station2); limit 99 stations
;     (TRI[i], I=1,NTR) (INT)   Caltech Notation
;                               10000*(I1)+100*(I2)+I3
; ** NOTE ERROR ALERT!
;   if one passes this 15 stations, but identifies one of the stations with
;   a number LARGER than 15, then there will be ERRORS if its merged with another
;   file.. So..
;        *** ALWAYS NUMBER THE STATIONS CONSECUTIVELY! ***
;               (i.e. no GAPS)
; **
;
;     PLUS the following lines repeated over the range of position angles
;
;     PA        (REAL)  Position Angle (Degrees).  This is normally set to
;                       Zero, but if non-zero then pretends the data was taken
;                       a later time T so that the sky has rotated the required
;                       PA.  This is a useful way to view data but only works
;                       if the same interferometric array is used.
;                       If one wanted to input EVERY speckle frame into this
;                       routine, a natural thing to do would be to
;                       have the PA advance, by using the Parallactic Angle found
;                       in the Keck Headers.
;     (AMP[I],        I=1,NBAS) (REAL)   Calibrated Visibility
;     (DELTAMP[I],    I=1,NBAS) (REAL)   Visibility Uncertainty
;     (PHASE[I],      I=1,NBAS) (REAL)   Dummy Phase.(deg)
;                       This will be passed with a 180 degree uncertainty, but
;                       appears here so that this routine will not have to do
;                       the Singular Value Decomposition, but rather lets the
;                       IDL routine come up with dummy phases which match the
;                       REAL closure phases.
;     (CLOS[I],       I=1,NTRI) (REAL)   Closure Phase (deg)
;     (DELTACLOS[I],  I=1,NTRI) (REAL)   Closure Phase Uncertainty (deg)
;
; 
;  However.. WHATS more important for YOU the USER is:
;    WHAT do I need to INPUT to this procedure.
;
;   If your data is from a non-redudant array, then this should be EASIER to DEAL
;   with.  If this is from a REDUNDANT array, then one must define MORE stations
;   than otherwise needed.
;
; keckmrg_input,x_coords,y_coords,Vis_array,vis_err_array,cp_array,cp_err_array,$
;    pa=pa,wave=wave,file=file,append=append,pref=pref,error_histogram=error_histogram
;
;  
;	x_coords(Ns)	Vector of X-coords for Ns stations. (meters)
;       y_coords(Ns) 	same for Y-coords (meters)
;	vis_array(ns,ns) This vector will contain all the visibility data 
;  		Since not all possible baselines may not exist in some cases,
;		a NEGATIVE value  for VIS_ERR_ARRAY 
;  		will indicate that this baseline dATA
;               does not EXIST.
;             NOTE: REQUIRES I1<I2 (I1,I2 are the staion numbers)
;	vis_err_arr(ns,ns) 	same as vis_array..
;       cp_array(ns,ns,ns)  Closure phase for Triangle ns,ns,ns.  Not all triangles
;       	 	are available so  (similar to vis_array), a triangle will
;    			be flagged as non-existent if the clp_err_array value is
; 			negative.
;		NOTE:  REQUIRES I1<I2<I3
;       clp_err_array(ns,ns,ns) same as clp_array
;
;	pa=position angle in degrees (default=0)
;       wave=wavelength (microns);  Default is 2.2
;       file=output file name.  Default: 'keckmrg_input.data'
;
;     OPTIONS:
;        /NOFIX    This suppresed the Closure Phase Calibration 
; 		   Routine by JDM.. This routine may introduce
;		   Artifacts into the data.
;
;	/PURGE	   This routine, when using NOFIX=0 will
;		   direct the FIX_CP routine to disregard
;		   closure triangles which don't make sense.
;		   That is, whose error bars on the reformulation are
;		   LARGER than the original one.. That is, either the
;		   uncertainty goes down or we eliminate it.
;	*  ** NOte: /PURGE is not implemented.
;
;	/APPEND    Will Append the new Angle to an existing file.
;
;       ERROR_HISTOGRAM  If greater than 2, this will plot clp errors and 
;                        prompt for cut. All all clp errors greater than
;                        the user input will be set to -1
;
; NOTE NOTE NOTE
;
;  I will write to the merge file EVERY baselines (i,j) such that
;  i<j.  Also every closing triangle (i,j,k) such that i<j<k.
;  However, if the baseline/triangle is invalid, the associated
;  error will always be -1.

pro keckmrg_input,x_coords,y_coords,Vis_array,vis_err_array,cp_array,cp_err_array,$
    pa=pa,wave=wave,file=file,nofix=nofix,purge=purge,append=append,pref=pref,    $
    error_histogram=error_histogram,ph_array=ph_array,$
    purge_baselines=purge_baselines

if (keyword_set(purge) eq 0) then purge=0
if (keyword_set(purge_baselines) eq 0) then purge_baselines=-1
if (keyword_set(pref) eq 0) then pref=0
if (keyword_set(error_histogram) eq 0) then error_histogram=0
if (keyword_set(x_coords) eq 0 ) then begin
print,'keckmrg_input,x_coords,y_coords,Vis_array,vis_err_array,cp_array,cp_err_array,
print,'    pa=pa,wave=wave,file=file,/nofix
return
endif
if (keyword_set(pa) eq 0) then pa=0.0
if (keyword_set(file) eq 0) then file='keckmrg_input.data'
if (keyword_set(wave) eq 0) then wave=2.2
if (keyword_set(nofix) eq 0) then nofix=0

nstats=long(n_elements(x_coords))
wav=wave*1e-6  ;convert to meters

; Now Create the cp and baseline index files.  Stations are labeled 1 and up.
baselines=lonarr(nstats,nstats)
cphases=lonarr(nstats,nstats,nstats)
for i=0l,nstats-2l do begin
  for j=i+1l,nstats-1l do begin
   baselines(i,j)=100l*(i+1l)+j+1l
  endfor
endfor
for i=0l,nstats-3l do begin
  for j=i+1l,nstats-2l do begin
    for k=j+1l,nstats-1l do begin
      cphases(i,j,k)=10000l*(i+1l)+100l*(j+1l)+k+1l
    endfor
  endfor
endfor

; Find non-zero baselines and closure triangles

baseindex1=where((vis_err_array ge 0.0) and (baselines gt 0l),nbas)
cpindex1=where((cp_err_array ge 0.0) and (cphases gt 0l),ntr)
baseindex=where((baselines gt 0l),nbas)
cpindex=where((cphases gt 0l),ntr)


  ; This ensures that I1<I2(<I3).

oldarray=cp_array
old_errors=cp_err_array
if (nofix eq 0) then begin
oldsigma=mean(cp_err_array(cpindex1))
print,'Beginning the Self-Consistency Calibration Procedure..'
fix_cp,cp_array,cp_err_array,change=change,purge=purge
print,'Fix CP Progress:  Iteration:',1,'  stdev(change):',stdev(change(cpindex1)),' degrees '
test_err=cp_err_array
newsigma=mean(cp_err_array(cpindex1))
tol=10.0
iter=1
while (tol gt 1.) do begin   ; Lets the stdev be less than 1 degree!
  iter=iter+1
fix_cp,cp_array,test_err,change=change,/keep
 tol=stdev(change(cpindex1))
print,'Fix CP Progress:  Iteration:',iter,'  stdev(change):',tol,' degrees '
endwhile
print,'  '
print,'Number of Iterations to FIX closure Phases: ',iter
print,'   Old <Closure Phase Error>: ',oldsigma
print,'   New <Closure Phase Error>: ',newsigma
print,'  '
endif
if (nofix gt 0) then begin   ; THIS MEANS THAT ONE IS SUPPRESSING THE
 			     ; FIX CP routine.
print,' '
print,' NOFIX IS ON:  WILL USE RAW Closure Phases.'
print,'   NOTE:  This means that the Phases passed can not be
print,'          made completely consistent with ALL the closure
print,'          phases.  Because the closure phases as they stand
print,'          are not internally consistent.'
print,'   The following histogram is the result of calculating a set of
print,'   phases from the existing closure phases and then recalculating
print,'   the closure phases.  This histogram shows the difference between
print,'   the original set of closure phases and the new set.'
; RESTORING.
cp_array=oldarray
cp_err_array=old_errors  

endif
choosecpagain:
if (error_histogram gt 1 ) then begin
   m=max(cp_err_array)
   plothist,cp_err_array(where(cp_err_array gt 0)),bin=.5, $
         xr=[0,m],title='Histogram of CLP errors'
   print,'Click on plot to define maximum error: points with large errors rejected'
   print,'Click to left of axis to exit'
   cursor,x,y & wait,.5
   print,x,y
  
   if(x gt 0.0 and x le max(cp_err_array(where(cp_err_array gt 0)))) then begin
     cp_err_array(where(cp_err_array gt x)) = -1.
   endif 
   plothist,cp_err_array(where(cp_err_array gt 0)),bin=.5, $
         xr=[0,m],title='Histogram of CLP errors'
endif


ph_array=float(calc_phases(cp_array,cp_err_array,pref=pref))

; JDM 29Nov98
if (purge_baselines eq 1) then begin
 print,'                   '
 print,'Now determining the minimal set of baseslines to invalidate
 print,'  In order to assure no bad triangles are included in the VLBI
 print,'  fits'       
 print,' '
return_bad=minset_badbaselines(ph_array,cp_err_array)
if (return_bad(0) ne -1) then begin
;  plothist,b_lengths(return_bad),bin=.4,tit="Histogram of Bad Baselines"
;  axaf=' ' & read,'Hit return to continue (or  r to choose new bad closure triangles: ',axaf 
;if (axaf eq 'r') then goto,choosecpagain 
vis_err_array(return_bad)=-1
endif
endif

   
;stop

; Check to make sure the inversion worked.
cp_test=fltarr(nstats,nstats,nstats)
for i=0,nstats-3 do begin
for j=i+1,nstats-2 do begin
 for k=j+1,nstats-1 do begin
  cp_test(i,j,k)=ph_array(i,j)+ph_array(j,k)-ph_array(i,k)
 endfor
endfor
endfor
if (n_elements(cpindex1) gt 1) then begin
print,'Results from PHASE calculation, derived from Closure Phases:'
print,'   Closure Phases: Standard Deviation of Differences'
print,'   ',stdev(angle_diff(cp_array(cpindex1),cp_test(cpindex1)))

q=uniq( angle_diff(cp_array(cpindex1),cp_test(cpindex1)) )
if (n_elements(q) gt 1) then begin
  plothist,angle_diff(cp_array(cpindex1),cp_test(cpindex1)),bin=.05
endif else print,'Exact Match.'
print,' '
endif
print,'Converting IDL Errors to VLBI Closure Phase Errors...'
cp_err_array=change_cp_errors(cp_err_array,bad_flag=45,bad_num=bad_num)
print,'  Removed ',bad_num,' Triangles with Raw Errors above 45 degrees.'
; Invalidates all Closure phases with idl error above (cp_err_cutoff) degrees.


print,'Now ready to WRITE'

if (keyword_set(append) eq 0) then $
openw,unit,file,/get,/swap_if_little_endian          else begin
  openu,unit,file,/get,/append,/swap_if_little_endian
  print," Appending to :",file
endelse

if (keyword_set(append) eq 0) then begin
printf,unit,nstats
printf,unit,wav
printf,unit,nbas
printf,unit,ntr
printf,unit,x_coords
printf,unit,y_coords
printf,unit,baselines(baseindex)
printf,unit,cphases(cpindex)
endif

printf,unit,pa
printf,unit,vis_array(baseindex)
printf,unit,vis_err_array(baseindex)
printf,unit,ph_array(baseindex)
printf,unit,cp_array(cpindex)
printf,unit,cp_err_array(cpindex)

close,unit
free_lun,unit

print,'File Written: ',file
end

