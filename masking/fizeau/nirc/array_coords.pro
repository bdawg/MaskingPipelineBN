; Feb 20, 1997
; jdm
;
; array_coords
;
; This function will return the (x,y,z) position in an array for a given
; set of indices.
;
;
;  indices needs to be 1-dimensionsal!
;
;  output of file is (NUMBER of INDICES, num of dim)
;   unless indices is a scalar, in which output is (num_dim)
;
;  example_array must be 
; 
;  function array_coords(indices,Example_array,help=help)
;
function array_coords,indices,example_array,help=help



if (keyword_set(help) eq 1) then begin
  print,' function array_coords(indices,Example_array,help=help)'
  return,0
endif

info1=size(example_array)
num=info1(0)

info2=size(indices)
if (info2(0) eq 0) then scalar=1 else scalar=info2(1)   ; must be either 0,1

init_string='result=lonarr('+string(scalar)
stars=' '
for i=0,num-1 do begin
  init_string=init_string+','+string(info1(i+1))
  stars=stars+',*'
endfor
  init_string=init_string+')'
;print,'Executing:  ',init_string
;q=execute(init_string)
result=lonarr(scalar,num)

for i=0,scalar-1 do begin
  index=indices(i)
  useful='result('+string(i)+stars+')= ['
  temp=fltarr(num)
  
  
     for dims=num-1,0,-1 do begin   ;Work Backwards.
        factor=1
        for kk=0,dims-1 do $
           factor=factor*info1(kk+1)
	temp(dims)=floor(index/factor)
        index=index-temp(dims)*factor
	result(i,dims)=temp(dims)
;        print,'i,factor,dims,temp(dims)
;	print,i,factor,dims,temp(dims)
     endfor
  useful=useful+string(temp(0))
    for lastloop=1,num-1 do useful=useful+','+string(temp(lastloop))
  useful=useful+']'
;  print,'Executing:  ',useful
;  q=execute(useful)
endfor
       
        


result=reform(result)
return,result
end

