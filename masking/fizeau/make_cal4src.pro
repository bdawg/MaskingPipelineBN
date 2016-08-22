;######## RELEASE CLIP HERE ######## 
;This script will make a cal4src matrix based on user interaction.

function make_cal4src, olog

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

print,  cal4src
return, cal4src

end 
