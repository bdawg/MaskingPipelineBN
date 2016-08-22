; JDM2005May30
; Document later.
; pass an array of times and a time internval.
; This routine will routine an index array that allow time averaging.
;
function timechunks,time,interval,nchunks=nchunks,counts=counts
; will return array of size (2, nchunks) = [start time, end time] for each
;				 	timechunk
; optional output variables:
;    nchunk = number of chunks found
;    counts = number of values in each chunk

time0=time(sort(time))
n=n_elements(time0)
counted=fltarr(n)

; START
tstart0=min(time0)
tend0 = tstart0+interval
in=where( time0 ge tstart0 and time0 lt tend0, ct)
counted(in)=1

tchunks=[ [tstart0,tend0] ]
counts=[ct]

; Now continue in loop
left=where(counted eq 0,ctleft)
while (ctleft gt 0) do begin
tstart0=min(time0(left))
tend0 = tstart0+interval
in=where( time0 ge tstart0 and time0 lt tend0, ct)
tchunks=[ [tchunks],[tstart0,tend0] ]
counted(in)=1
counts=[counts,ct]

left=where(counted eq 0,ctleft)
endwhile
info=size(tchunks,/dim)
if n_elements(info) eq 1 then nchunks=1 else nchunks=info(1)

return,tchunks
end



