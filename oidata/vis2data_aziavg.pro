; 2005May30 JDM
; pass oivis2data array and length. This program will
; bin data in chunks of baseline radius and average all vis2data quantities,
; returning a new one.
; 
; This really wants an OITABLE (from read_oidata) NOT a idl vis2data
; structure (from extract_vis2data)
;
; 
; Will respect different 
; 	target_ids
;	insname 
;
; Required helper: timechunks.pro
function vis2data_azIAvg, oivis2,interval,$ ; interval in meters
	extra1=ex1,out_extra1=out_extra1

; will do proper AVERAGE of quantity ex1. e.g. Hour Angle
ucall=oivis2.ucoord
vcall=oivis2.vcoord
ri2at,ucall,vcall,aaall,ttall


targets=oivis2(uniq(oivis2.target_id,sort(oivis2.target_id))).target_id
nt=n_elements(targets)

insnames = oivis2(uniq(oivis2.insname,sort(oivis2.insname))).insname
ni=n_elements(insnames)

 for t0=0,nt-1 do begin
  for i0=0,ni-1 do begin


 in_original = where (   $
	      oivis2.target_id eq targets(t0) and $
              oivis2.insname eq insnames(i0) , ct0)

    if (ct0 gt 0) then begin

counted=fltarr(ct0)
dummy=where(counted eq 0,num_left)
  
while (num_left gt 0) do begin
  in0=where (  counted ne 1,ct)
  in=in_original(in0)
  uc=ucall(in)
  vc=vcall(in)
  ri2at,uc,vc,aa,tt
  instart=(where(aa eq min(aa)))(0)
;for vis2data, reflect,(u,v)?. ; do not do this for t3data/visdata
;  diff2a=(uc-uc(instart))^2 + (vc-vc(instart))^2
;  diff2b=(uc+uc(instart))^2 + (vc+vc(instart))^2
;  diff2 = diff2a<diff2b
diff2=abs(aa-aa(instart))
  ingood=where(diff2 lt interval,ct)
  in=in_original(in0(ingood))
  counted(in0(ingood))=1 ;sorry for this confusing method...

  if ct eq 1 then begin
   new_tab=oivis2(in)
   if keyword_set(ex1) then new_ex1=ex1(in)
  endif ; ct=1

  if (ct gt 1) then begin
     subset=oivis2(in)
     new_tab=subset(0)
     new_tab.time = median(subset.time)
     new_tab.int_time = max(subset.time) - min(subset.time)
     new_tab.mjd = median(subset.mjd)
 ; make sure you get the MEAN baseline and set PA to 0.0
     ri2at,subset.ucoord,subset.vcoord,subsetaa,subsetpa
     new_tab.ucoord=mean(subsetaa)
     new_tab.vcoord=0

   if keyword_set(ex1) then new_ex1=mean(ex1(in))
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
     if keyword_set(ex1) then out_ex1=new_ex1
  endif else begin
     new_oivis2=concat_oitable(new_oivis2,new_tab)
      if keyword_set(ex1) then out_ex1=[out_ex1,new_ex1]
  endelse

;counted     
dummy=where(counted eq 0,num_left)
endwhile

endif ; at least one chunk
endfor
endfor
if (keyword_set(ex1)) then out_extra1=out_ex1
return,new_oivis2

end

