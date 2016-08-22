; JDM
;  ri2at.pro
;
;   REAL/Imaginary to Amplitude/Phase
;

pro ri2at,re,im,am,ph,help=help

if (keyword_set(help) eq 1) then begin
  print,' REAL/Imaginary to Amplitude/Phase (phase in Degrees)
 print,'pro ri2at,re,im,am,ph
 return
endif
pi=3.14159265358979323846264338327950288419716939937d
;stop
ph=re
ph(*)=0.0
am=abs(complex(re,im))
in=where(re ne 0,ct)
if (ct ne 0) then $
  ph(in)=atan(im(in)/re(in))*180.d /pi
in=where(re lt 0,count)
if (count gt 0) then $
  ph(in)=ph(in) + 180.d
ph=(ph+360.0000d) mod 360.d

index=where(re eq 0 and im lt 0.0,ct) 
if (ct gt 0) then begin
  ph(index)=270.d
  am(index)=abs(im(index))
endif
index=where(re eq 0 and im gt 0.0,ct)
if (ct gt 0) then begin
  ph(index)=90d
  am(index)=abs(im(index))
endif


end


