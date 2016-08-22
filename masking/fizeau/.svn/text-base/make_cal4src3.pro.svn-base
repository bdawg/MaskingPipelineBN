;######## RELEASE CLIP HERE ######## 
;This script will make a cal4src matrix based on user interaction.

function make_cal4src3, cubeinfo, quad=quad

restore,cubeinfo

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (olog.instrument[0] eq 'CONICA') then begin

nt = n_elements(olog.tsize) 

;; Indices corresponding to unique values of the FULL olog.tsize array with
;; the values arranged in increasing order:
s = sort(olog.tsize) 
u =  uniq(olog.tsize, s)

;; Now to accommodate for cases when either every source_name is
;; recorded AND also those cases where only the source_name's
;; corresponding to each unique tsize value are recordedd:
ixs=uniq(olog.tsize,sort(olog.tsize))
utsize=olog.tsize[ixs]
usnames=olog.source_name[ixs]
ucnames=olog.cube_fname[ixs]
incr=sort(utsize)

print,  'Choose blocks to analyse together (enter -1 when finished...)' 
for i = 0, n_elements(u)-1 do begin 
  print,strtrim(i, 2)+'  '+olog.cube_fname[min(where(olog.tsize eq olog.tsize[u[i]]))]+'   '+strcompress(usnames[incr[i]],/remove_all)
endfor 
temp=1
tunum=intarr(1) 
tnums=intarr(1)
while temp ge 0 do begin
  read,  temp
  tunum=[tunum,temp]
  if temp ge 0 then tnums=[tnums,where(olog.tsize eq olog.tsize[u[temp]])]
endwhile
tunum=tunum[1:*]
tnums=tnums[1:*]
tunum=tunum[where(tunum ge 0)]
src=u[tunum]
 
cunum=intarr(1)
dummy=0
for i=0,n_elements(u)-1 do begin
  for j=0,n_elements(tunum)-1 do begin
    if olog.tsize[u[i]] eq olog.tsize[u[tunum[j]]] then dummy=1
  endfor
  if dummy ne 1 then cunum=[cunum,i]
  dummy=0
endfor
cunum=cunum[1:*]
cals=u[cunum]

while (n_elements(cals) gt 1) do begin
  print,  'Choose calibrators to reject (-1 for none).'
  for i = 0, n_elements(cals)-1 do begin 
    print,  strtrim(i,2)+'  '+ucnames[where(utsize eq olog.tsize[cals[i]])]+'   '+usnames[where(utsize eq olog.tsize[cals[i]])]
  endfor
  read,  crej
  if (crej eq -1) then break
  cals =  cals(where(cals ne cals[crej]))
endwhile


cal4src =  intarr(nt, nt)
for i = 0, n_elements(cals)-1 do begin
 cnums =  where(olog.tsize eq olog.tsize[cals[i]])
 for j = 0, n_elements(tnums)-1 do cal4src[cnums,tnums[j]] = 1
endfor

while (n_elements(tnums) gt 1) do begin
 rej=''
 print,  'Num  Name    Frame   Tsize '
 print,  '--------------------------------------------'
 for i = 0, n_elements(tnums)-1 do begin 
;  print,  tnums[i], olog.cube_fname[tnums[i], 1],olog.tsize[tnums[i]],olog.utc[tnums[i]], format = '(I-4, A-8, I-4,A-12)'
  print,  tnums[i], usnames[where(utsize eq olog.tsize[tnums[i]])], '   ', olog.frames[tnums[i]],'   ', olog.tsize[tnums[i]], format = '(I-4, A-8, A-3, I-5, A-3, I-4)'
 endfor
 print,  'Choose a target cube to reject or'
 print,  'enter first and last cubes for a contiguous block (i.e. first-last)'
 print,  '(-1 for none)'
 read,  rej
 if (rej[0] eq '-1') then break
 pos=strpos(rej,'-')
 if pos lt 0 then begin 
   cal4src[*, rej] = 0
   tnums =  tnums[where(tnums ne rej)]
 endif else begin
   first=fix(strmid(rej,0,pos))
   last=fix(strmid(rej,pos+1,strlen(rej)-pos))
   rej=indgen(last-first+1)+first
   for i=0,n_elements(rej)-1 do begin
     cal4src[*,rej[i]]=0
     tnums=tnums[where(tnums ne rej[i])]
   endfor
  endelse
endwhile
 
;; If requested, remove from the analysis any frames not taken with
;; the source in the specified chip quadrant:

if keyword_set(quad) then begin

;  w = where(framestats.quad ne quad)
;  for i=0,n_elements(w)-1 do begin
;    cal4src[w[i],*]=0
;    cal4src[*,w[i]]=0
;  endfor

  for i=0,nt-1 do begin
    check = 0
    for j=0,n_elements(quad)-1 do begin
      if framestats.quad[i] eq quad[j] then check=1
    endfor
    if check eq 0 then begin
      cal4src[i,*]=0
      cal4src[*,i]=0
    endif
  endfor

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endif else begin

if keyword_set(quad) then begin
  print,'Quadrant selection only available for CONICA currently'
  print,'To continue,using all quadrants, enter .c'
  stop
endif

nt = n_elements(olog.source_name)
s = sort(olog.source_name)
u =  uniq(olog.source_name, s)
print,  'Choose the target to analyse...' 
for i = 0, n_elements(u)-1 do begin 
 print,  strtrim(i, 2), ' ', olog.source_name(u[i])
endfor 
read,  tunum ;Unique target number
cals =  where(u ne u[tunum])
tnums =  where(olog.source_name eq olog.source_name[u[tunum]])
;For now, we will choose calibrators by name. Trimming the array further comes later.
while (n_elements(cals) gt 1) do begin
 print,  'Choose a calibrator to reject (-1 for none).'
 for i = 0, n_elements(cals)-1 do begin 
  print,  strtrim(cals[i], 2), ' ', olog.source_name(u[cals[i]])
 endfor 
 read,  crej
 if (crej eq -1) then break
 cals =  cals(where(cals ne crej))
endwhile

cal4src =  intarr(nt, nt)
for i = 0, n_elements(cals)-1 do begin
 cnums =  where(olog.source_name eq olog.source_name[u[cals[i]]])
 for j = 0, n_elements(tnums)-1 do cal4src[cnums, tnums[j]] = 1
endfor

while (n_elements(tnums) gt 1) do begin
 print,  'Num Cube            UT'
 print,  '--------------------------------------------'
 for i = 0, n_elements(tnums)-1 do begin 
  print,  tnums[i], olog.cube_fname[tnums[i], 1],olog.utc[tnums[i]], format = '(I-4, A-16, A-12)'
 endfor
 print,  'Choose a target cube to reject (-1 for none)'
 read,  rej
 if (rej eq -1) then break
 cal4src[*, rej] = 0
 tnums =  tnums(where(tnums ne rej))
endwhile

endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return, cal4src

end 








