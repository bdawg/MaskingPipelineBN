;######## RELEASE CLIP HERE ######## 
; Program to flag bad data frames for Conica
;
; Version 0.0 From cube_nirc2  PGT 25Oct05

function flagbad_conica,fstats,cube,discard_sigma,tot_bad

  spk_sz=size(fstats)
  nspeck=spk_sz(2) 

  tot_bad=0

; %%%%% MANUAL bad frame flagging
  if(discard_sigma[0] ne -1) then begin
  if(discard_sigma[0] eq 0) then begin
    n_bad=0
    bad_frames=-1
    good_frames=indgen(nspeck)
    print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print,'You have selected the option to manually discard frames'
    print,'To accept a frame, click right of the y axis'
    print,'To reject a frame, click left  of the y axis'
    print,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    for i=0,nspeck-1 do begin
      image_cont,cube(*,*,i),/nocont,/asp,tit=string(i)+' click >0 accept <0 reject' 
      cursor,x,y & wait,.4
      if(x lt 0) then begin
         if(n_bad eq 0) then in=[i] else in=[in,i]
         n_bad=n_bad+1
      endif
    endfor
    if ((n_bad) ge nspeck) then begin
      print,'ERROR:  ALL Speckle Frames are flagged bad'
      stop
    endif
    if (n_bad gt 0) then begin
      bad_frames=in
      good_frames(bad_frames)=-1
      good_index=good_frames(where(good_frames ge 0))
      print,''
      print,'### Script will REMOVE ### ',n_bad,' Bad Frames from your arrays.'
      print,''
      tot_bad=n_bad
      ;frames=frames(*,*,good_index)
    endif

  endif else begin  ;%%%%% AUTOMATIC bad frame flagging
    if (nspeck le 3) then begin
      print,'You only have ',nspeck,' frames ... are you sure you want to do statistics?' 
      stop
    endif
;   Now Flag all speckle sets which are (discard_sigma) standard deviations
;   from the mean (in ANY statistic (except SKY LEVEL))
    resistant_mean,fstats(0,*),2.0,mn0 & sd0=robust_sigma(fstats(0,*))
    resistant_mean,fstats(1,*),2.0,mn1 & sd1=robust_sigma(fstats(1,*))
    resistant_mean,fstats(2,*),2.0,mn2 & sd2=robust_sigma(fstats(2,*))
    sd2 = sd2 >  0.03*mn2 ;;Sometimes resistant_mean has a silly answer for small number of frames
    resistant_mean,fstats(3,*),2.0,mn3 & sd3=robust_sigma(fstats(3,*))
    sd3 = sd3 >  0.06*mn3 ;;Sometimes resistant_mean has a silly answer for small number of frames
    md3 =  median(fstats[3, *])

; If discard_sigma is a 5 element vector, implement the AO-dropout code
    if ( (size(discard_sigma))[1] eq 5) then aotop=discard_sigma[4] else aotop=-1
 

; Find outliers - given the discard_sigma is set positive
    if(discard_sigma[0] gt 0) then begin
      in0=where(  ( fstats(0,*) gt (mn0+discard_sigma(0)*sd0) ) or  $
	        ( fstats(0,*) lt (mn0-discard_sigma(0)*sd0) ) ,ct0) 
    endif else begin
      in0=-1 & ct0=0
    endelse
    
    if(discard_sigma[1] gt 0) then begin
    in1=where(  ( fstats(1,*) gt (mn1+discard_sigma(1)*sd1) ) or  $
                ( fstats(1,*) lt (mn1-discard_sigma(1)*sd1) ) ,ct1)
    endif else begin
      in1=-1 & ct1=0
    endelse
    
    if(discard_sigma[2] gt 0) then begin
    in2=where(  ( fstats(2,*) gt (mn2+discard_sigma(2)*sd2) ) or  $
                ( fstats(2,*) lt (mn2-discard_sigma(2)*sd2) ) ,ct2)
    endif else begin
      in2=-1 & ct2=0
    endelse
    
    if(discard_sigma[3] gt 0 and aotop le 0) then begin
    in3=where(  ( fstats(3,*) gt (mn3+discard_sigma(3)*sd3) ) or  $
                ( fstats(3,*) lt (mn3-discard_sigma(3)*sd3) ) ,ct3)
    endif else begin
      in3=-1 & ct3=0
    endelse

    ; Remove MJI cut on fractional peak flux. More sophisticated one below.
    ;in4 = where(fstats[3, *] lt discard_sigma[4]*md3)

; -- This is AO Dropout code which removes frames based on stats of aotop fraction --
; %%%% Discard cut based on peak counts to get rid of low-strehl dropouts
; %%%% This is done by computing statistics of only a certain fraction
; %%%%%  discard_sigma[4] %  of the highest pk pixel images.
    if(aotop gt 0) then begin
      st=sort(fstats(3,*))
      pksort=fstats(3,st)
      nsamp=nint(nspeck*discard_sigma[4])   ; Grab this fraction
      hi_pk_samp=pksort[nspeck-nsamp:nspeck-1] ; of the highest pk flx
      ;resistant_mean,hi_pk_samp,2.0,samp_mn
      samp_mn=median(hi_pk_samp)
      samp_sd=robust_sigma(hi_pk_samp)
      lowercut=samp_mn - discard_sigma[3] * samp_sd
      in4=where(fstats(3,*) lt lowercut, ct4)
      
  ;     mn=samp_mn & sg=samp_sd*2.0
  ;     plot,fstats[3,*]
  ;     oplot,[0,100],[mn,mn]
  ;     oplot,[0,100],[mn-sg,mn-sg],line=1
  ;     oplot,[0,100],[mn+sg,mn+sg],line=1
  ;     oplot,(indgen(100))(in4),(fstats[3,*])(in4),psym=4
  ;     stop
    endif else begin
      in4 = -1
      ct4 = 0 
    endelse

    in=[in0,in1,in2,in3, in4]

;  %%%%% Flag bad frames
    n_bad=0
    bad_frames=-1
    good_frames=indgen(nspeck)
    if ( (ct0+ct1+ct2+ct3+ct4) gt 0) then begin
      in=in(where(in ne -1))
      in=in(sort(in))
      bad_frames=[in(uniq(in))]
      n_bad=n_elements(bad_frames)
      if ((n_bad) ge nspeck-2) then begin
        print,'ERROR:  ALL Speckle Frames are flagged bad'
        stop
      endif
      good_frames(bad_frames)=-1
      good_index=good_frames(where(good_frames ge 0))
      print,''
      print,'### Script will REMOVE ### ',n_bad,' Bad Frames from your arrays.'
      print,''
      tot_bad=n_bad
      ;cube=cube(*,*,good_index)
    endif
   endelse
   endif else begin  ;%%%%% Case: don't even check anymore
    n_bad=0
    bad_frames=-1
    good_frames=indgen(nspeck)
   endelse
; %%%%% Finsihed Discard Frames

return,good_frames

end
