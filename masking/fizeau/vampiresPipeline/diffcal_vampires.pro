; Reads bs_....idlvar files and creates differentially calibrated observables.
;
; This early version just does simple division, treating each pair of
; channels as orthogonal. Future version will use the proper measured
; Mueller matrices for each channel.
;
; BN Aug2013

; Version BN Nov2013 - Implement Closure Phase

; -------------------------- Input Files -----------------------
;set_plot,'ps'
;device,filename='vhvv_vega18h775_single_hwp.ps',ysize=28,xsize=20,bits=8,/color,yoffset=0,xoffset=0


pro diffcal_vampires, fileprefix, startnum, nfiles, cubeinfofile, csvoutfile, fluxvec=fluxvec,$
                      yrange=yrange, nbootstraps=nbootstraps, docals=docals, saveeps=saveeps, $
                      output_special=output_special



prefix=fileprefix; 'VCalTests_ann_775_turb0_pQWPs0_switching_seq_'
;output_special=''

if ~keyword_set(nbootstraps) then nbootstraps=100
if ~keyword_set(docals) then DoCals = ['0', '1a', '2a', '3a', '4a', '5a', '6a']
if ~keyword_set(saveeps) then saveeps=1
if ~keyword_set(output_special) then output_special = ''

; HWP=0 degrees files:
;hwp0FileNums=[24,28]

; HWP=22.5 degrees files:
;hwp225FileNums=[25,29]

; HWP=45 degrees files:
;hwp45FileNums=[26,30]

; HWP=67.5 degrees files:
;hwp675FileNums=[27,31]

;makefilenums,0,64,hwp0FileNums,hwp225FileNums,hwp45FileNums,hwp675FileNums
makefilenums,startnum,nfiles,hwp0FileNums,hwp225FileNums,hwp45FileNums,hwp675FileNums

; cubeinfo file:
;;;;cubeinfoFile='cubeinfoAug2014.idlvar'


; ---------------------------- Options --------------------------
longestBL = -1;2e6 ;Ignore baselines longer than this (-1 for none)

pairwise = 0 ;Set to 1 to do divisions before averaging.
             ;Perhaps useful for highly correlated pairs (eg Wollaston channels)

;calType = '6a'  ; 0 = Triple Calibration
                ; 1a = Double cal, Woll. + LCVR, HWP0
                ; 1b = Double cal, Woll. + LCVR, HWP45
                ; 2a = Double cal, Woll + HWP, LCVR1
                ; 2b = Double cal, Woll + HWP, LCVR2
                ; 3a = Double cal, LCVR + HWP, Woll.1
                ; 3b = Double cal, LCVR + HWP, Woll.2
                ; 4a = Single cal, Wollaston (LCVR1, HWP0)
                ; 4b = Single cal, Wollaston (LCVR1, HWP45)
                ; 4c = Single cal, Wollaston (LCVR2, HWP0)
                ; 4d = Single cal, Wollaston (LCVR2, HWP45)
                ; 5a = Single cal, LCVR (Woll1, HWP0)
                ; 6a = Single cal, HWP (LCVR1, Woll1)


if keyword_set(yrange) eq 0 then yrange=[0.9,1.1]


;Qonly = 0    ; Set to 1 to only do Stokes Q (ie HWP=0,45)

bootstraps = nbootstraps    ; No. of bootstrap iterations to find variances
                     ; Set to 1 for no error calcs.

clip0 = 1            ; Set to 1 to clip negative V2s to zero.

subtmean = 0         ; Set to 1 to subtract mean VH/VV
;----------------------------------------------------------------

if keyword_set(fluxvec) eq 0 then fluxvec=fltarr(4)

ncals=n_elements(docals)
AllSDs=fltarr(ncals, 2) ;For Q and U

for ii = 0,ncals-1 do begin
calType=DoCals[ii]






; Read in all data:
readdiffbs, prefix, hwp0FileNums, hwp225FileNums, hwp45FileNums, hwp675FileNums, $
            h0_v2s_all, h0_bss_all, h225_v2s_all, h225_bss_all, $
            h45_v2s_all, h45_bss_all, h675_v2s_all, h675_bss_all, $
            u_coords, v_coords

restore,cubeinfoFile

lamstring=olog.filter
lamstring=strsplit(lamstring[0],'-',/extract) ;Assumes centre wavelength before '-'
reads,lamstring[0],lambda
lambda = lambda*1e-9 ;Put in m

if clip0 eq 1 then begin
    h0_v2s_all=h0_v2s_all>0
    h225_v2s_all=h225_v2s_all>0
    h45_v2s_all=h45_v2s_all>0
    h675_v2s_all= h675_v2s_all>0
endif

if longestBL[0] gt 0 then begin
    blengths = sqrt(u_coords^2+v_coords^2)
    gd=where(blengths LT longestBL)
    h0_v2s_all=h0_v2s_all[gd,*,*,*]
    h225_v2s_all=h225_v2s_all[gd,*,*,*]
    h45_v2s_all=h45_v2s_all[gd,*,*,*]
    h675_v2s_all=h675_v2s_all[gd,*,*,*]
    u_coords=u_coords[gd]
    v_coords=v_coords[gd]
endif


; Set up bootstrapping
    nbls=(size(h0_v2s_all))(1)
    ncps=(size(h0_bss_all))(1)
if bootstraps gt 1 then begin
    h0_v2s_all_orig=h0_v2s_all
    h225_v2s_all_orig=h225_v2s_all
    h45_v2s_all_orig=h45_v2s_all
    h675_v2s_all_orig=h675_v2s_all
    vhvv_bootstraps=dblarr(nbls,bootstraps)
    vhvvU_bootstraps=dblarr(nbls,bootstraps)
    nf=(size(h0_v2s_all))(2)
    h0_bss_all_orig=h0_bss_all
    h225_bss_all_orig=h225_bss_all
    h45_bss_all_orig=h45_bss_all
    h675_bss_all_orig=h675_bss_all
    cp_bootstraps=dblarr(ncps,bootstraps)
    cpU_bootstraps=dblarr(ncps,bootstraps)
endif

for bootstrap_it = 0,bootstraps-1 do begin

if bootstraps gt 1 then begin
    resamp_inds=floor(randomu(seed,nf)*nf)
    h0_v2s_all=h0_v2s_all_orig[*,resamp_inds,*,*]
    h45_v2s_all=h45_v2s_all_orig[*,resamp_inds,*,*]
    h0_bss_all=h0_bss_all_orig[*,resamp_inds,*,*]
    h45_bss_all=h45_bss_all_orig[*,resamp_inds,*,*]
endif 


; Triple calibration:
if pairwise eq 1 then begin
    ; Here, do divisions *before* averaging.
    ; Perhaps useful for highly correlated pairs (eg Wollaston channels)
    print,'Not yet implemented, sorry.'
    stop
endif else begin
    ; Average then divide - the normal way.

    ; Do Stokes Q

    ; 'H' and 'V' comments are just to improve readbility, they're
    ; not actually H and V polarisations.
    ; Convention used is that Woll. ch1 by itself would be H.
    ; NB LCVR Voltage 1 is 1/2 wave retardance, Voltage 2 is ~0 retardance.
    h0_l1_w1=total(h0_v2s_all[*,*,0,0],2) / (size(h0_v2s_all))(2)    ; V
    h0_l1_w2=total(h0_v2s_all[*,*,1,0],2) / (size(h0_v2s_all))(2)    ; H
    h0_l2_w1=total(h0_v2s_all[*,*,0,1],2) / (size(h0_v2s_all))(2)    ; H
    h0_l2_w2=total(h0_v2s_all[*,*,1,1],2) / (size(h0_v2s_all))(2)    ; V
    h45_l1_w1=total(h45_v2s_all[*,*,0,0],2) / (size(h45_v2s_all))(2)    ; H
    h45_l1_w2=total(h45_v2s_all[*,*,1,0],2) / (size(h45_v2s_all))(2)    ; V
    h45_l2_w1=total(h45_v2s_all[*,*,0,1],2) / (size(h45_v2s_all))(2)    ; V
    h45_l2_w2=total(h45_v2s_all[*,*,1,1],2) / (size(h45_v2s_all))(2)    ; H
   
    h0_lcvr1 = sqrt( h0_l1_w1 / h0_l1_w2 )    ; (V/H)^1/2
    h0_lcvr2 = sqrt( h0_l2_w1 / h0_l2_w2 )    ; (H/V)^1/2
    h45_lcvr1 = sqrt( h45_l1_w1 / h45_l1_w2 )    ; (H/V)^1/2
    h45_lcvr2 = sqrt( h45_l2_w1 / h45_l2_w2 )    ; (V/H)^1/2

    h0 = sqrt( h0_lcvr1 / h0_lcvr2 ) ; (V/H)^1/2
    h45 = sqrt( h45_lcvr1 / h45_lcvr2 ) ; (H/V)^1/2
;;h45 = sqrt( h45_lcvr2 / h45_lcvr1 ) ; (H/V)^1/2

    vhvv = sqrt( h45 / h0 )    ; (H/V)^1/2
                               ; NB this ^1/2 is from inherited convention,
                               ; wherein observable is V_H/V_V, not V^2_H/V^2_V.


    ; Now do Closure Phases
    cp_h0_l1_w1=atan(total(h0_bss_all[*,*,0,0],2),/phase)    ; V
    cp_h0_l1_w2=atan(total(h0_bss_all[*,*,1,0],2),/phase)    ; H
    cp_h0_l2_w1=atan(total(h0_bss_all[*,*,0,1],2),/phase)    ; H
    cp_h0_l2_w2=atan(total(h0_bss_all[*,*,1,1],2),/phase)    ; V
    cp_h45_l1_w1=atan(total(h45_bss_all[*,*,0,0],2),/phase)    ; H
    cp_h45_l1_w2=atan(total(h45_bss_all[*,*,1,0],2),/phase)    ; V
    cp_h45_l2_w1=atan(total(h45_bss_all[*,*,0,1],2),/phase)    ; V
    cp_h45_l2_w2=atan(total(h45_bss_all[*,*,1,1],2),/phase)    ; H

    cp_h0_lcvr1 = ( cp_h0_l1_w1 - cp_h0_l1_w2 )    ; V-H
    cp_h0_lcvr2 = ( cp_h0_l2_w1 - cp_h0_l2_w2 )    ; H-V
    cp_h45_lcvr1 = ( cp_h45_l1_w1 - cp_h45_l1_w2 )    ; H-V
    cp_h45_lcvr2 = ( cp_h45_l2_w1 - cp_h45_l2_w2 )    ; V-H

    cp_h0 = ( cp_h0_lcvr1 - cp_h0_lcvr2 )/2 ; [(V-H)-(H-V)]/2 = V-H
    cp_h45 = ( cp_h45_lcvr1 - cp_h45_lcvr2 )/2 ; [(H-V)-(V-H)]/2 = H-V

    cp = (cp_h0 - cp_h45)/2 ;[(V-H)-(H-V)]/2 = V-H


    ; Do other calibration types
    case caltype of
        '0' : ;Do nothing
        '1a': begin
              vhvv = h0
              cp = cp_h0
              end
        '1b': begin
              vhvv = h45
              cp = cp_h45
              end
        '2a': begin
              vhvv = sqrt( h0_lcvr1 / h45_lcvr1 )
              cp = (cp_h0_lcvr1 - cp_h45_lcvr1)/2
              end
        '2b': begin
              vhvv = sqrt( h0_lcvr2 / h45_lcvr2 )
              cp = (cp_h0_lcvr2 - cp_h45_lcvr2)/2
              end
        '3a': begin
              top=sqrt( h0_l1_w1 / h0_l2_w1 )
              bot=sqrt( h45_l1_w1 / h45_l2_w1 )
              vhvv=sqrt(top/bot)
              top=(h0_l1_w1 - h0_l2_w1)
              bot=(h45_l1_w1 - h45_l2_w1)
              cp=fltarr(n_elements(cp_h0)) & print,'CP Not Yet Implemented for this cal mode'
              end
        '3b': begin
              top=sqrt( h0_l1_w2 / h0_l2_w2 )
              bot=sqrt( h45_l1_w2 / h45_l2_w2 )
              vhvv=sqrt(top/bot)
              top=(h0_l1_w2 - h0_l2_w2)
              bot=(h45_l1_w2 - h45_l2_w2)
              cp=fltarr(n_elements(cp_h0)) & print,'CP Not Yet Implemented for this cal mode'
              end
        '4a': begin
              vhvv = h0_lcvr1
              cp = cp_h0_lcvr1
              end
        '4b': begin
              vhvv = h45_lcvr1
              cp = cp_h45_lcvr1
              end
        '4c': begin
              vhvv = h0_lcvr2
              cp = cp_h0_lcvr2
              end
        '4d': begin
              vhvv = h45_lcvr2
              cp = cp_h45_lcvr2
              end
        '5a': begin
              vhvv = sqrt( h0_l1_w1 / h0_l2_w1 )
              cp = cp_h0_l1_w1 - cp_h0_l2_w1
              end
        '6a': begin
              vhvv = sqrt( h0_l1_w1 / h45_l1_w1 )
              cp = cp_h0_l1_w1 - cp_h45_l1_w1
              end
    endcase


   if bootstraps gt 1 then begin
       vhvv_bootstraps[*,bootstrap_it]=vhvv
       cp_bootstraps[*,bootstrap_it]=cp
   endif

    ; #### TO DO - Closure phases!
endelse

endfor

if bootstraps gt 1 then begin
    vhvv = dblarr(nbls)
    vhvverr = dblarr(nbls)
    cp = dblarr(ncps)
    cperr = dblarr(ncps)
    for kk = 0,nbls-1 do begin
        vhvv[kk]=mean(vhvv_bootstraps[kk,*])
        vhvverr[kk]=stddev(vhvv_bootstraps[kk,*])
    endfor
    for kk = 0,ncps-1 do begin
        cp[kk]=mean(cp_bootstraps[kk,*])
        cperr[kk]=stddev(cp_bootstraps[kk,*])
    endfor
endif else begin
    ;vhvv=vhvv_bootstraps[*,0]
    vhvverr=dblarr(nbls)
    cperr=dblarr(ncps)
endelse

;stop
; Plot results
blengths = sqrt(u_coords^2+v_coords^2)*lambda
bazims = atan(v_coords/u_coords)
!p.multi=[0,1,4]
!p.charsize=3.0;1.6
loadct,39,/silent

if subtmean eq 1 then vhvv-=mean(vhvv)-1

ploterr,blengths,vhvv,vhvverr,title=['Standard deviation: '+strn(stddev(vhvv))],yr=yrange,ytit='Polarised Visibility Ratio',xtit='Baseline length at pupil (m)',psym=4
oplot,[0,1e8],[1,1],linestyle=1

plot,bazims,vhvv,xtitle='Baseline azimuth angle (rads)',ytitle='Polarised Visibility Ratio',$
  psym=4,/nodata,yrange=yrange
for j=0,n_elements(bazims)-1 do oploterr,[bazims[j],bazims[j]],[vhvv[j],vhvv[j]],$
  [vhvverr[j],vhvverr[j]], color=(blengths[j]/max(blengths))*250,psym=4
oplot,[-2,2],[1,1],linestyle=1




; ------------------------------- Repeat for HWP225-675 -----------------------------
for bootstrap_it = 0,bootstraps-1 do begin

if bootstraps gt 1 then begin
    resamp_inds=floor(randomu(seed,nf)*nf)
    h225_v2s_all=h225_v2s_all_orig[*,resamp_inds,*,*]
    h675_v2s_all=h675_v2s_all_orig[*,resamp_inds,*,*]
endif 


; Triple calibration:
if pairwise eq 1 then begin
    ; Here, do divisions *before* averaging.
    ; Perhaps useful for highly correlated pairs (eg Wollaston channels)
    print,'Not yet implemented, sorry.'
    stop
endif else begin
    ; Average then divide - the normal way.

    ; Do Stokes Q

    ; 'H' and 'V' comments are just to improve readbility, they're
    ; not actually H and V polarisations.
    ; Convention used is that Woll. ch1 by itself would be H.
    ; NB LCVR Voltage 1 is 1/2 wave retardance, Voltage 2 is ~0 retardance.
    h225_l1_w1=total(h225_v2s_all[*,*,0,0],2) / (size(h225_v2s_all))(2)    ; V
    h225_l1_w2=total(h225_v2s_all[*,*,1,0],2) / (size(h225_v2s_all))(2)    ; H
    h225_l2_w1=total(h225_v2s_all[*,*,0,1],2) / (size(h225_v2s_all))(2)    ; H
    h225_l2_w2=total(h225_v2s_all[*,*,1,1],2) / (size(h225_v2s_all))(2)    ; V
    h675_l1_w1=total(h675_v2s_all[*,*,0,0],2) / (size(h675_v2s_all))(2)    ; H
    h675_l1_w2=total(h675_v2s_all[*,*,1,0],2) / (size(h675_v2s_all))(2)    ; V
    h675_l2_w1=total(h675_v2s_all[*,*,0,1],2) / (size(h675_v2s_all))(2)    ; V
    h675_l2_w2=total(h675_v2s_all[*,*,1,1],2) / (size(h675_v2s_all))(2)    ; H
   
    h225_lcvr1 = sqrt( h225_l1_w1 / h225_l1_w2 )    ; (V/H)^1/2
    h225_lcvr2 = sqrt( h225_l2_w1 / h225_l2_w2 )    ; (H/V)^1/2
    h675_lcvr1 = sqrt( h675_l1_w1 / h675_l1_w2 )    ; (H/V)^1/2
    h675_lcvr2 = sqrt( h675_l2_w1 / h675_l2_w2 )    ; (V/H)^1/2

    h225 = sqrt( h225_lcvr1 / h225_lcvr2 ) ; (V/H)^1/2
    h675 = sqrt( h675_lcvr1 / h675_lcvr2 ) ; (H/V)^1/2

    vhvvU = sqrt( h675 / h225 )    ; (H/V)^1/2
                               ; NB this ^1/2 is from inherited convention,
                               ; wherein observable is V_H/V_V, not V^2_H/V^2_V.

    ; Do other calibration types
    case caltype of
        '0' : ;Do nothing
        '1a': vhvvU = h225
        '1b': vhvvU = h675
        '2a': vhvvU = sqrt( h225_lcvr1 / h675_lcvr1 )
        '2b': vhvv = sqrt( h225_lcvr2 / h675_lcvr2 )
        '3a': begin
              top=sqrt( h225_l1_w1 / h225_l2_w1 )
              bot=sqrt( h675_l1_w1 / h675_l2_w1 )
              vhvvU=sqrt(top/bot)
              end
        '3b': begin
              top=sqrt( h225_l1_w2 / h225_l2_w2 )
              bot=sqrt( h675_l1_w2 / h675_l2_w2 )
              vhvvU=sqrt(top/bot)
              end
        '4a': vhvvU = h225_lcvr1
        '4b': vhvvU = h675_lcvr1
        '4c': vhvvU = h675_lcvr2
        '4d': vhvvU = h675_lcvr2
        '5a': vhvvU = sqrt( h225_l1_w1 / h225_l2_w1 )
        '6a': vhvvU = sqrt( h225_l1_w1 / h675_l1_w1 )
    endcase


   if bootstraps gt 1 then begin
       vhvvU_bootstraps[*,bootstrap_it]=vhvvU
   endif

    ; #### TO DO - Closure phases!
endelse

endfor

if bootstraps gt 1 then begin
    vhvvU = dblarr(nbls)
    vhvvUerr = dblarr(nbls)
    for kk = 0,nbls-1 do begin
        vhvvU[kk]=mean(vhvvU_bootstraps[kk,*])
        vhvvUerr[kk]=stddev(vhvvU_bootstraps[kk,*])
    endfor
endif else begin
    ;vhvv=vhvv_bootstraps[*,0]
    vhvvUerr=dblarr(nbls)
endelse


; Plot results

if subtmean eq 1 then vhvvU-=mean(vhvvU)-1

ploterr,blengths,vhvvU,vhvvUerr,title=['Standard deviation: '+strn(stddev(vhvvU))],yr=yrange,ytit='Polarised Visibility Ratio',xtit='Baseline length at pupil (m)',psym=4
oplot,[0,1e8],[1,1],linestyle=1

plot,bazims,vhvvU,xtitle='Baseline azimuth angle (rads)',ytitle='Polarised Visibility Ratio',$
  psym=4,/nodata,yrange=yrange
for j=0,n_elements(bazims)-1 do oploterr,[bazims[j],bazims[j]],[vhvvU[j],vhvvU[j]],$
  [vhvvUerr[j],vhvvUerr[j]], color=(blengths[j]/max(blengths))*250,psym=4
oplot,[-2,2],[1,1],linestyle=1


;device,/close
;set_plot,'x'


save,blengths,bazims,vhvv,vhvverr,vhvvu,vhvvuerr,file=['diffdata_'+prefix+'_'+output_special+caltype+'.idlvar']


;;;;;;;;;;;;; Output plots
if saveeps eq 1 then begin
   plotfile=['diffdata_'+prefix+'_'+output_special+caltype+'.eps']
   comments=[prefix+'_'+output_special+caltype+'.idlvar']

   set_plot,'ps'
   device,filename=plotfile,ysize=18,xsize=12,bits_per_pixel=8,/color,yoffset=0,xoffset=0,/encapsulated

   !p.multi=[0,1,4]
   !p.charsize=1.6
   loadct,39,/silent
;stop

   ploterr,blengths,vhvv,vhvverr,title=['Standard deviation: '+strn(stddev(vhvv))],yr=yrange,ytit='Polarised Visibility Ratio',xtit='Baseline length at pupil (m)',psym=4
   oplot,[0,1e8],[1,1],linestyle=1

   plot,bazims,vhvv,xtitle='Baseline azimuth angle (rads)',ytitle='Polarised Visibility Ratio',$
        psym=4,/nodata,yrange=yrange,title=comments
   for j=0,n_elements(bazims)-1 do oploterr,[bazims[j],bazims[j]],[vhvv[j],vhvv[j]],$
                                            [vhvverr[j],vhvverr[j]], color=(blengths[j]/max(blengths))*250,psym=4
   oplot,[-2,2],[1,1],linestyle=1

   ploterr,blengths,vhvvU,vhvvUerr,title=['Standard deviation: '+strn(stddev(vhvvU))],yr=yrange,ytit='Polarised Visibility Ratio',xtit='Baseline length at pupil (m)',psym=4
   oplot,[0,1e8],[1,1],linestyle=1

   plot,bazims,vhvvU,xtitle='Baseline azimuth angle (rads)',ytitle='Polarised Visibility Ratio',$
        psym=4,/nodata,yrange=yrange,title=comments
   for j=0,n_elements(bazims)-1 do oploterr,[bazims[j],bazims[j]],[vhvvU[j],vhvvU[j]],$
                                            [vhvvUerr[j],vhvvUerr[j]], color=(blengths[j]/max(blengths))*250,psym=4
   oplot,[-2,2],[1,1],linestyle=1

   device,/close
   set_plot,'x'
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plotlinkstring = ['=hyperlink("./'+plotfile+'")']
openw,1,csvoutfile,/append
outstr=[fileprefix+'_'+output_special+' , '+caltype+' , '+strn(stddev(vhvv))+' , '+strn(stddev(vhvvU))+' , '+strn(mean(vhvverr))+' , '+strn(mean(vhvvUerr))+' , '+strn(fluxvec[0])+' , '+strn(fluxvec[1])+' , '+strn(fluxvec[2])+' , '+strn(fluxvec[3])+' , ' + plotlinkstring]
printf,1,outstr
close,1



endfor

skip:
return

end
;; wset,1
;; !p.multi=0
;; plothist,cp/!pi*180,/auto
;; print,stddev(cp/!pi*180)
;; wset,0
;; end
