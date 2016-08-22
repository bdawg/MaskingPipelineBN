; plota_nirc2.pro
; Generates plots of raw data products for NIRC2
; data analysis suite. Based on get_pow_hc.pro
; from the original Nirc data software suite.  
; (Version: most recent comment Sep01)
;
; Original Version                              PGT     08Dec03
; Nirc2 Version                                 PGT     28Feb04
;
; To be done : add plots of raw speckle data frames.

pro plota_nirc2,stats,date,headinfo,num,history,comments,			$
	a_destripe_flag=a_destripe_flag,saturation_flag=saturation_flag,        $
	hardcopy=hardcopy,bad_frames=bad_frames,identifier=identifier


if (keyword_set(stats) eq 0) then begin
  print,'pro  plota_nirc,stats,date,headinfo,num,history,comments...'
  return
endif
if (keyword_set(identifier) eq 0) then identifier='nosave'

; %%%%% Information from headinfo data structure ...
name=headinfo.source_name[0]
color=headinfo.filter[0]
mask='' 	;ignore for now
dimx=headinfo.nax1[0]
dimy=headinfo.nax2[0]
pa=headinfo.pa[0]
	
printloop=0
cleanplot,/silent

secondloop4print:

if (keyword_set(saturation_flag) eq 0) then saturation_flag=0
if (keyword_set(a_destripe_flag) eq 0) then a_destripe_flag=0
if ((size(bad_frames))[0] eq 0) then num_bad=0 else num_bad=n_elements(bad_frames) 

!p.multi=[0,2,2]
plot,stats(0,*),stats(1,*),tit='Center of Speckle Cloud',$
   xtit="X Axis (Azimuth)",$
   ytit="Y Axis (Elevation)",psym=3,xr=[0,dimx-1],yr=[0,dimy-1],xstyle=1,ystyle=1
if (num_bad gt 0) then $
  oplot,[stats(0,bad_frames)],[stats(1,bad_frames)],psym=4
legend,[date,name+'('+num+')',color+'/'+mask,history,comments,$
  'Saturate Flag='+strtrim(string(nint(saturation_flag)),2)+$
'; A_Destripe Flag='+strtrim(string(nint(a_destripe_flag)),2)+$
'; Bad Frames='+strtrim(string(nint(num_bad)),2)+$
'; PA='+strtrim(string(nint(pa)),2)+' Deg'],$
  line=[-1,-1,-1,-1,-1,-1],box=0,pos=[.05,.95],spacing=.09,charsize=.7,$
  /norm

center_x=median(stats(0,*))
center_y=median(stats(1,*))
angs=7.*findgen(10000)/9999.
sigx=robust_sigma(stats(0,*))
sigy=robust_sigma(stats(1,*))
plots,center_x+sigx*cos(angs),center_y+sigy*sin(angs),psym=3

ii=findgen(  (size(stats))(2))
plot,stats(2,*),title='Integrated Flux',xtitle='Speckle Frame #'
if (num_bad gt 0) then $
  oplot,[ii(bad_frames)],[stats(2,bad_frames)],psym=4
sd=robust_sigma(stats(2,*))
resistant_mean,stats(2,*),2.0,mn
legend,[strtrim(string(mn),2)+' +/- '+strtrim(string(sd),2)],pos=[.52,.65],$
  box=0,/norm,line=-1
oplot,[0,1000],[mn,mn],noclip=0
oplot,[0,1000],[mn+sd,mn+sd],noclip=0,line=2
oplot,[0,1000],[mn-sd,mn-sd],noclip=0,line=2
 
plot,stats(3,*),title='Peak Flux',xtitle='Speckle Frame #'
if (num_bad gt 0) then $
  oplot,[ii(bad_frames)],[stats(3,bad_frames)],psym=4
sd=robust_sigma(stats(3,*))
resistant_mean,stats(3,*),2.0,mn
legend,[strtrim(string(mn),2)+' +/- '+strtrim(string(sd),2)],pos=[.02,.15],box=0,$
 /norm,line=-1
 
oplot,[0,1000],[mn,mn],noclip=0
oplot,[0,1000],[mn+sd,mn+sd],noclip=0,line=2
oplot,[0,1000],[mn-sd,mn-sd],noclip=0,line=2


plot,stats(4,*),title='Sky Background Offset',xtitle='Speckle Frame #'
if (num_bad gt 0) then $
  oplot,[ii(bad_frames)],[stats(4,bad_frames)],psym=4
sd=robust_sigma(stats(4,*))
resistant_mean,stats(4,*),2.0,mn

legend,[strtrim(string(mn),2)+' +/- '+strtrim(string(sd),2)],line=-1,pos=[.52,.15],box=0,$
  /norm
oplot,[0,1000],[mn,mn],noclip=0
oplot,[0,1000],[mn+sd,mn+sd],noclip=0,line=2
oplot,[0,1000],[mn-sd,mn-sd],noclip=0,line=2

if (identifier ne 'nosave' or keyword_set(hardcopy) eq 1) then begin
   set_plot,'ps'
   device,/landscape
   printloop=printloop+1
endif

if(printloop eq 1) then goto,secondloop4print

if (keyword_set(hardcopy) ne 0) then begin
  device,/close
  spawn,'lpr idl.ps'
endif
if (identifier ne 'nosave') then begin
  device,/close
  spawn,'\mv idl.ps '+identifier
endif
set_plot,'x'

cleanplot,/silent
end



