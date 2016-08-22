; Program to flag bad data frames for Nirc2
;
; Version 0.0 From cube_nirc2  PGT 25Oct05
;;If the first element of discard_sigma is -1... don't do
;;anything. Before it wouldn't do anything if this was negative.
;;                             MJI 26March2010

function flagbad_nirc2,frames,fstats,discard_sigma,filenumbers,  tot_bad

  spk_sz=size(frames)
  if(spk_sz[0] gt 2) then nspeck=spk_sz(3) else nspeck=1

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
              filenum=["n"+string(filenumbers[i],format="(I4.4)")+".fits"]
              image_cont,frames(*,*,i),/nocont,/asp,tit=filenum+' click >0 accept <0 reject' 
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
              print,'### YOU NEED TO REMOVE ### ',n_bad,' Bad Frames from your arrays.'
              print,' Run script over with src files:',filenumbers[good_index] 
              print,''
              tot_bad=tot_bad+1
                                ;frames=frames(*,*,good_index)
          endif
          
      endif else begin          ;%%%%% AUTOMATIC bad frame flagging
          if (nspeck le 3) then begin
              print,'You only have ',nspeck,' frames ... are you sure you want to do statistics?' 
              stop
          endif
;   Now Flag all speckle sets which are (discard_sigma) standard deviations
;   from the mean (in ANY statistic (except SKY LEVEL))
          resistant_mean,fstats(0,*),2.0,mn0 & sd0=robust_sigma_mike(fstats(0,*))
          resistant_mean,fstats(1,*),2.0,mn1 & sd1=robust_sigma_mike(fstats(1,*))
          resistant_mean,fstats(2,*),2.0,mn2 & sd2=robust_sigma_mike(fstats(2,*))
          sd2 = sd2 >  0.03*mn2 ;;Sometimes resistant_mean has a silly answer for small number of frames
          resistant_mean,fstats(3,*),2.0,mn3 & sd3=robust_sigma_mike(fstats(3,*))
          sd3 = sd3 >  0.06*mn3 ;;Sometimes resistant_mean has a silly answer for small number of frames
          md3 =  median(fstats[3, *])
          
;;Only count image motions beyond 2 pixels (MJI)
          in0=-1 & in1=-1 & in2=-1 & in3=-1 & in4=-1
          ct0=0 & ct1=0 & ct2=0 & ct3=0 & ct4=0
          if (discard_sigma[0] gt 0) then in0=where(  ( fstats(0,*)-mn0 gt (discard_sigma(0)*sd0 > 2) ) or  $
                      ( mn0-fstats(0,*) gt (discard_sigma(0)*sd0 > 2) ) ,ct0)
          if (discard_sigma[1] gt 0) then in1=where(  ( fstats(1,*)-mn1 gt (discard_sigma(1)*sd1 > 2) ) or  $
                      ( mn1-fstats(1,*) gt (discard_sigma(1)*sd1 > 2) ) ,ct1)
          if (discard_sigma[2] gt 0) then in2=where(  ( fstats(2,*) gt (mn2+discard_sigma(2)*sd2) ) or  $
                      ( fstats(2,*) lt (mn2-discard_sigma(2)*sd2) ) ,ct2)
          if (discard_sigma[3] gt 0) then in3=where(  ( fstats(3,*) gt (mn3+discard_sigma(3)*sd3) ) or  $
                      ( fstats(3,*) lt (mn3-discard_sigma(3)*sd3) ) ,ct3)
          if (discard_sigma[4] gt 0) then in4 = where(fstats[3, *] lt discard_sigma[4]*md3, ct4)
          in=[in0,in1,in2,in3,in4]
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
              print,'### YOU NEED TO REMOVE ### ',n_bad,' Bad Frames from your arrays.'
              print,' Run script over with src files:',filenumbers[good_index] 
              print,''
              tot_bad=tot_bad+1
                                ;cube=cube(*,*,good_index)
          endif
      endelse
  endif else begin              ;%%%%% Case: don't even check anymore
      n_bad=0
      bad_frames=-1
      good_frames=indgen(nspeck)
  endelse
; %%%%% Finsihed Discard Frames

return,good_frames

end
