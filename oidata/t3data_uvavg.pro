; 2005May30 JDM
; pass oit3data array and a time scale. This program will
; bin data in chunks of time and average all t3data quantities,
; returning a new one.
; **
; * I know that technically one should average the t3amplitude like
; a vector (bispectrum). however, since we have errors on all the
; amp and phase separatley, I will average separately in case
; t3amp is not a good observable [technically it should be flagged as such,
; but...]
; 
; This really wants an OITABLE (from read_oidata) NOT a idl t3data
; structure (from extract_t3data)
; 
; Will respect different 
; 	target_ids
;	insname 
;       sta_indexes
;       arrname
;
; Required helper: timechunks.pro
function t3data_uvavg, oit3,interval,$ ; interval in meters
        extra1=ex1,out_extra1=out_extra1

; will do proper AVERAGE of quantity ex1. e.g. Hour Angle

ucall1=oit3.u1coord
vcall1=oit3.v1coord
ri2at,ucall1,vcall1,aaall1,ttall1
ucall2=oit3.u2coord
vcall2=oit3.v2coord
ri2at,ucall2,vcall2,aaall2,ttall2


targets=oit3(uniq(oit3.target_id,sort(oit3.target_id))).target_id
nt=n_elements(targets)

insnames = oit3(uniq(oit3.insname,sort(oit3.insname))).insname
ni=n_elements(insnames)

arrnames = oit3(uniq(oit3.arrname,sort(oit3.arrname))).arrname
na=n_elements(arrnames)

sta0 = oit3(uniq(oit3.sta_index[0],sort(oit3.sta_index[0]))).sta_index[0]
ns0=n_elements(sta0)

sta1 = oit3(uniq(oit3.sta_index[1],sort(oit3.sta_index[1]))).sta_index[1]
ns1=n_elements(sta1)

sta2 = oit3(uniq(oit3.sta_index[2],sort(oit3.sta_index[2]))).sta_index[2]
ns2=n_elements(sta2)

 for t0=0,nt-1 do begin
  for i0=0,ni-1 do begin
   print,t0+1,nt,i0+1,ni
   for a0=0,na-1 do begin
    for s0=0,ns0-1 do begin
     print,a0+1,na,s0+1,ns0
     for s1=0,ns1-1 do begin
      for s2=0,ns2-1 do begin

 in_original = where (   $
              oit3.target_id eq targets(t0) and $
              oit3.insname eq insnames(i0) and $
              oit3.arrname eq arrnames(a0) and $
              oit3.sta_index[0] eq sta0(s0) and $
              oit3.sta_index[1] eq sta1(s1) and $
              oit3.sta_index[2] eq sta2(s2), ct0)

    if (ct0 gt 0) then begin

counted=fltarr(ct0)
dummy=where(counted eq 0,num_left)

while (num_left gt 0) do begin
  in0=where (  counted ne 1,ct)
  in=in_original(in0)
  uc1=ucall1(in)
  vc1=vcall1(in)
  uc2=ucall2(in)
  vc2=vcall2(in)
  ri2at,uc1,vc1,aa1,tt1
  ri2at,uc2,vc2,aa2,tt2

  instart=(where(aa1 eq min(aa1)))(0)
;for vis2data, reflect,(u,v)?. ; do not do this for t3data/visdata
;  diff2a=(uc-uc(instart))^2 + (vc-vc(instart))^2
;  diff2b=(uc+uc(instart))^2 + (vc+vc(instart))^2
;  diff2 = diff2a<diff2b
diff1=(uc1-uc1(instart))^2 + (vc1-vc1(instart))^2
diff2=(uc2-uc2(instart))^2 + (uc2-uc2(instart))^2

  ingood=where(diff2 lt interval*interval and diff1 lt interval*interval,ct)
  in=in_original(in0(ingood))
  counted(in0(ingood))=1 ;sorry for this confusing method...


  if ct eq 1 then begin
   new_tab=oit3(in)
   if (keyword_set(ex1)) then new_ex1=ex1(in)
  endif ; ct=1


  if (ct gt 1) then begin
     subset=oit3(in)
     new_tab=subset(0)
     new_tab.time = median(subset.time)
     new_tab.int_time = max(subset.time) - min(subset.time)
     new_tab.mjd = median(subset.mjd)
     new_tab.u1coord=mean(subset.u1coord)
     new_tab.v1coord=mean(subset.v1coord)
     new_tab.u2coord=mean(subset.u2coord)
     new_tab.v2coord=mean(subset.v2coord)
     if (keyword_set(ex1)) then new_ex1=mean(ex1(in))
 
;---------------------------------------------------------------
; First do for t3amp
;---------------------------------------------------------------
; now average data and do error analysis
     nwave=new_tab.nwave
     vis2=fltarr(ct,nwave) ; re-using name -- really t3amp
     err= fltarr(ct,nwave)
     flag = bytarr(ct,nwave)
      ; get data
      for v=0,ct-1 do begin
       tempv=*subset(v).t3amp
       vis2(v,*)=tempv ; re-using names..... not really vis2 but t3amp
       tempe=*subset(v).t3amperr
       err(v,*)=tempe
       tempf=*subset(v).flag
       flag(v,*)=tempf
      endfor
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
    new_tab.t3amp=ptr_new(tempv)
    new_tab.t3amperr =ptr_new(tempe)
    new_tab.flag = ptr_new(tempf)
;---------------------------------------------------------------
; First do for t3phi
;---------------------------------------------------------------
; now average data and do error analysis
     nwave=new_tab.nwave
     vis2=fltarr(ct,nwave) ; re-using name -- really t3amp
     err= fltarr(ct,nwave)
     flag = bytarr(ct,nwave)
      ; get data
      for v=0,ct-1 do begin
       tempv=*subset(v).t3phi
       vis2(v,*)=tempv ; re-using names..... not really vis2 but t3amp
       tempe=*subset(v).t3phierr
       err(v,*)=tempe
       tempf=*subset(v).flag
       flag(v,*)=tempf
      endfor
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
           result=wtd_mean_angle(v2(goodindex),e2(goodindex),werr)
           ;vector averaging
           tempv(v)=result
           tempe(v)=werr
           tempf(v)=(byte('F'))(0)
	;*** Sanity Check ***
	; sometimes scatter between points dominates over the
	;errors bars, rendering the wtd mean inacccurate.
        ; Thus this routine will check to make sure the chi2/dof is one
	; after averaging.  If yes, the nothing. if chi2/dof > 1 then
	; this will multiple errors by factor to force chi2/dof=1.
         adiff=angle_diff(v2,tempv)
         chi2=total(   (adiff/e2(goodindex))^2)/$
		(goodct-1)
	  if (chi2 gt 1.25) then begin
           factor=sqrt(chi2)
           tempe(v)=tempe(v)*factor
                   endif 
        endif ; goodct =1
     endfor
    new_tab.t3phi=ptr_new(tempv)
    new_tab.t3phierr =ptr_new(tempe)
    new_tab.flag = ptr_new(tempf)
;---------------------------------------------------------------
endif ;ct>1

  if (n_elements(new_oit3) eq 0) then begin
     new_oit3=new_tab
     if (keyword_set(new_ex1)) then out_ex1=new_ex1
  endif else begin
     new_oit3=concat_oitable(new_oit3,new_tab)
     if (keyword_set(new_ex1)) then out_ex1=[out_ex1,new_ex1]
  endelse


;counted
dummy=where(counted eq 0,num_left)
endwhile
     
endif; nchunks>0
endfor
endfor
endfor
endfor
endfor
endfor
if (keyword_set(out_ex1)) then out_extra1=out_ex1
return,new_oit3

end

