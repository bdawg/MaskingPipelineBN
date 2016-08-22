; february 20, 1997
;  
;  28Oct97 JDM
;    antisymmetrize phase array
;    Added /check feature
;
; JDM
;
; calc_phases
;
;  This routine is based on the method found in PHASES.FOR in the
;  caltech vlbi analysis package.  I have don't have any nice features, like
;  filling in phases to fit a model.  Until the number of phase reach a critical
;  mass, the phases will be filled in to equal zero.  also, there is presently
;  NO Order of Preference in filling in baselines.  It would make sense to choose
;  short baselines to be zero, especially for Diameters.. but.. not 
;  implemented
;
;  One will want to have a MODEL some points, and
;  also it would be useful to have this Randomly go through the phases
;  each time so that the iterative procedure doesn't get locked into particular
;  phases due a particular Fixed set of Baselines which are GIVEN from the model
;  each time.
;
;  The syntax is similar to that for svd_phases , but doesn't not require you
;  to pass a baseline array, since this will return a completely filled in
;  phase array, even for a sparsely passed triangle array.
;  You must pass this an array of Closure Phases (nb,nb,nb).
;  Also you must pass an array of closure phase errors, although
;  it is only use to look for values of -1, which means to disregard the
;  closure phase.
;
;  This routine will only use those which have a non-negative
;  value in the cp_err_array and
;  also only those which satisfy that (I1<i2(<I3))
;

;
;  phases=calc_phases(closure_phases)
;
function calc_phases,cp_array,cp_err_array,random=random,$
   model_phases=model_phases,preferences=preferences,check=check

info=size(cp_array)
nstats=long(info(1))
if (keyword_set(model_phases) eq 0) then model_phases=fltarr(nstats,nstats)
; we must prepare a preference matrix.  (0 for most preferred to use model,
;		-1 if not part of model)

;
if (keyword_set(preferences) eq 0) then begin
  preferences=replicate(-1,nstats,nstats)
  number=nstats*(nstats-1)/2.  ; This is an inelegant way to get 
			       ; a randomn preference matrix.

  rnd=randomu(seed,number)
  order=sort(rnd)
  c=0
  for i=0,nstats-2 do begin   
   for j=i+1,nstats-1 do begin 
    preferences(i,j)=order(c)
    c=c+1
   endfor
  endfor
endif
num_set=0
phases=fltarr(nstats,nstats)
known=fltarr(nstats,nstats)  ;To keep track of when  a baseline phase is determined.
unusable=fltarr(nstats,nstats,nstats)
index=where(cp_err_array eq -1,ct)
if (ct ne 0) then unusable(index)=1.0  ; unusable triangles 

model_phase=0.0   ; The phase of the unknown

; First loop through the Baselines.. i<j

list=0

for dummy_i=0,nstats-2 do begin
 for dummy_j=dummy_i+1,nstats-1 do begin
  ; Which baseline is this for? (must be legitimate one (i<j)
   pref_index=where(preferences eq list,howmany)
   if (howmany eq 0) then print,'ERROR: All preferences must be assigned'
			
   sothen=array_coords(pref_index,preferences)
   if (howmany eq 0) then print,'you are screwed.  All preferences must be given ( at the moment)'
   i=sothen(0)
   j=sothen(1)
   list=list+1
   if (known(i,j) eq 1) then goto,nextphase   ;Skip Triangles if this Baseline is known
  
   phases(i,j)=model_phases(i,j)
   known(i,j)=1
   num_set=num_set+1

; Print Status Report
 
num_known=0
num_total=0
for iii=0,nstats-2 do begin
 for jjj=iii+1,nstats-1 do begin
  num_total=num_total+1
  if (known(iii,jjj) eq 1) then num_known=num_known+1
  endfor
endfor
print,num_known*100./num_total,'% of Phases are now Fixed.'

moretriangles:
tricount=0
; Next we must loop through all the triangles and
; fill the phases if 2 or more are known.
   for k=0,nstats-3 do begin
    for l=k+1,nstats-2 do begin
     for m=l+1,nstats-1 do begin
  if (unusable(k,l,m) eq 0) then begin    ; Only go on if USABLE
      knowns=[known(k,l),known(l,m),known(k,m)]
      if (total(knowns) eq 3) then begin
		unusable(k,l,m)=1
;		print,'This point should never be reached'
		goto,nexttriangle  ; This should never happen
      endif

      if (total(knowns) eq 2) then begin  ; We can FIX one..
	tricount=tricount+1
	unusable(k,l,m)=1.0
        bindex=(where(knowns eq 0))(0)  ;Figure out which phase is unknown
        if (bindex eq 0) then begin 
 		phases(k,l)=cp_array(k,l,m)-phases(l,m)+phases(k,m)
 		known(k,l)=1.0
  	endif
	if (bindex eq 1) then begin 
		phases(l,m)=cp_array(k,l,m)-phases(k,l)+phases(k,m)
		known(l,m)=1.0	
	endif
	if (bindex eq 2) then begin 
		phases(k,m)=-1.0*cp_array(k,l,m)+phases(k,l)+phases(l,m)
                known(k,m)=1.0
        endif
      endif   ; This IF referred to whether or not there were 2 unknowns!
endif  ; This if referred to whether or not the triangle was usuable.
nexttriangle:
      endfor
    endfor
   endfor  ;We've looped through all the Triangles
if (tricount gt 0) then goto, moretriangles ; Since we had luck last round, do it
					    ; again. If no Triangles were found
					    ; to be completed, then set the
					    ;present phase to its model values
					    ; and continue on.


nextphase:
  endfor
endfor


; All phases should be known now..
; Check to see!
;

for i=0,nstats-2 do begin
 for j=i+1,nstats-1 do begin
  if (known(i,j) eq 0) then begin
    print,'Lordie! There is a big bug in here! Better stop now.'
	stop
   endif
  endfor
endfor


phases=mod360(phases)
; Anti-Symmetrize the Phase Matrix for Plotting Purposes
for i=0,nstats-2 do begin
 for j=i+1,nstats-1 do begin
   phases(j,i)=-1*phases(i,j)
endfor
endfor
print,'Number of Phases SET: ',num_set

if (keyword_set(check) eq 1) then begin
; Check to make sure the inversion worked.
 
cp_test=fltarr(nstats,nstats,nstats)
for i=0,nstats-3 do begin
for j=i+1,nstats-2 do begin
 for k=j+1,nstats-1 do begin
  cp_test(i,j,k)=phases(i,j)+phases(j,k)-phases(i,k)
 endfor
endfor
endfor
test_cp,qq,qqq,num=nstats,in=cpindex1
cpindex1=cpindex1(where(cp_err_array(cpindex1) ne -1))
print,'Results from PHASE calculation, derived from Closure Phases:'
print,'   Closure Phases: Standard Deviation of Differences'
print,'   ',stdev(angle_diff(cp_array(cpindex1),cp_test(cpindex1)))


plothist,angle_diff(cp_array(cpindex1),cp_test(cpindex1)),bin=.05

endif


return,phases
end






