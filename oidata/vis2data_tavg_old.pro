; 2005May30 JDM
; pass oivis2data array and a time scale. This program will
; bin data in chunks of time and average all vis2data quantities,
; returning a new one.
; 
; This really wants an OITABLE (from read_oidata) NOT a idl vis2data
; structure (from extract_vis2data)
;
; 
; Will respect different 
; 	target_ids
;	insname 
;       sta_indexes
;       arrname
;
; Required helper: timechunks.pro
function vis2data_tavg, oivis2,interval ; tint in minutes.


times=oivis2.mjd*24.*60. ; in minutes
tchunks=timechunks(times,interval,nch=nch)

targets=oivis2(uniq(oivis2.target_id,sort(oivis2.target_id))).target_id
nt=n_elements(targets)

insnames = oivis2(uniq(oivis2.insname,sort(oivis2.insname))).insname
ni=n_elements(insnames)

arrnames = oivis2(uniq(oivis2.arrname,sort(oivis2.arrname))).arrname
na=n_elements(arrnames)

sta0 = oivis2(uniq(oivis2.sta_index[0],sort(oivis2.sta_index[0]))).sta_index[0]
ns0=n_elements(sta0)

sta1 = oivis2(uniq(oivis2.sta_index[1],sort(oivis2.sta_index[1]))).sta_index[1]
ns1=n_elements(sta1)

for chunks=0,nch-1 do begin
 for t0=0,nt-1 do begin
  for i0=0,ni-1 do begin
   for a0=0,na-1 do begin
    for s0=0,ns0-1 do begin
     for s1=0,ns1-1 do begin

 in = where (   $
     times ge tchunks(0,chunks) and times lt tchunks(1,chunks) and $
              oivis2.insname eq insnames(i0) and $
              oivis2.arrname eq arrnames(a0) and $
              oivis2.sta_index[0] eq sta0(s0) and $
              oivis2.sta_index[1] eq sta1(s1), ct)
  if ct eq 1 then begin
   new_tab=oivis2(in)
  endif ; ct=1

  if (ct gt 1) then begin
     subset=oivis2(in)
     new_tab=subset(0)
     new_tab.time = mean(subset.time)
     new_tab.int_time = max(subset.time) - min(subset.time)
     new_tab.mjd = mean(subset.mjd)
     new_tab.ucoord=mean(subset.ucoord)
     new_tab.vcoord=mean(subset.vcoord)
; now average data and do error analysis
     nwave=new_tab.nwave
     vis2=fltarr(ct,nwave)
     err= fltarr(ct,nwave)
     flag = bytarr(ct,nwave)
      ; get data
      for v=0,ct-1 do begin
       tempv=*subset(v).vis2data
       vis2(v,*)=tempv
       tempe=*subset(v).vis2err
       err(v,*)=tempe
       tempf=*subset(v).flag
       flag(v,*)=tempf
      endfor
;print,*oivis2(0).vis2data
;stop
    ; Now do each element of data vector..
     for v=0,nwave -1 do begin
	
        v2=reform(vis2(*,v))
        e2=reform(err(*,v))
        f2=reform(flag(*,v))
        goodindex=where(e2 gt 0 and f2 eq (byte('F'))(0),goodct)
        if goodct eq 0 then begin
           tempv(v)=0.0
           tempe(v)=-1
           tempf(v)=(byte('T'))(0)
	endif
	if goodct eq 1 then begin
           tempv(v)=v2(goodindex)
           tempe(v)=e2(goodindex)
           tempf(v)=f2(goodindex)
        endif
        if (goodct gt 1) then begin
           tempv(v)=wtd_mean(v2(goodindex),e2(goodindex),newerr)
           tempe(v)=newerr
           tempf(v)=(byte('F'))(0)
	;*** Sanity Check ***
	; sometimes scatter between points dominates over the
	;errors bars, rendering the wtd mean inacccurate.
        ; Thus this routine will check to make sure the chi2/dof is one
	; after averaging.  If yes, the nothing. if chi2/dof > 1 then
	; this will multiple errors by factor to force chi2/dof=1.
         chi2=total(  ((v2(goodindex)-tempv(v))/(e2(goodindex)) )^2)/$
		(goodct-1)
	  if (chi2 gt 1.25) then begin
           factor=sqrt(chi2)
           tempv(v)=wtd_mean(v2(goodindex),e2(goodindex)*factor,newerr)
           tempe(v)=newerr
          endif 
        endif ; goodct =1
     endfor
    new_tab.vis2data=ptr_new(tempv)
    new_tab.vis2err =ptr_new(tempe)
    new_tab.flag = ptr_new(tempf)

endif ;ct>1

  if (n_elements(new_oivis2) eq 0) then begin
     new_oivis2=new_tab
  endif else new_oivis2=concat_oitable(new_oivis2,new_tab)


     
endfor
endfor
endfor
endfor
endfor
endfor
return,new_oivis2

end

