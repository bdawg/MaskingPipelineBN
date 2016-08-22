;+
; Pass aline of text and it will seperate into a vector of strings, each one
; was separated by white space
; if NOEQUAL is turned on then it first sets all equal signs to blanks
; if CLEANUP is on, it turns all ,=$ to blanks.
;
;function string_vector,st,count,cleanup=cleanup,del=del
;-
function string_vector,st,count,cleanup=cleanup,del=del

if (keyword_set(del) eq 0) then del=' '
a=st
if (keyword_set(cleanup) eq 1) then begin 
for i=0,strlen(a)-1 do begin
if (strmid(a,i,1) eq "=") then strput,a,' ',i
if (strmid(a,i,1) eq "$") then strput,a,' ',i
if (strmid(a,i,1) eq ",") then strput,a,' ',i

endfor
endif


a=strcompress(strtrim(a,2))
place=fltarr(100)
place(0)=-1
i=0
repeat begin
i=i+1
p=strpos(a,del,place(i-1)+1)
place(i)=p
endrep until (p eq -1)
place(i)=strlen(a)
place=place(0:i)

vect=strarr(i)
for j=0,i-1 do begin
  vect(j)=strmid(a,place(j)+1,(place(j+1)-place(j)-1))
endfor
count=i

return,vect
end

