; This .pro is based on the older cube_vampires.script

; This takes FITS files from VAMPIRES and produces cubes
; (and associated cubeinfo file) suitable for further processing
; with the IDL aperture-masking pipeline.
;
; To maintain compatability, the two LCVR states are treated as two
; separate cubes - i.e. each raw FITS file now produces 4 cubes. This is
; recorded in the lcvrstate variable, as a value of 1 or 2.
;
; BS files will have suffixes '_1' or '_2' for up or down respectively,
; AND '_A' or '_B' for LCVR state 1 or state 2 respectively.

; Next step e.g.:
; calc_bispect,'cubeinfoTestingAug2013.idlvar' , root_dir='~/code_svn/masking/',plotall=2,ddir='./'

pro cube_vampires
  
; All settings that used to be in cube_vampires.script are now in this
; common block:
common procsettings


if (cal_start[0] ge 0) then begin 
  frames = [src_start+indgen(nsrc), cal_start+indgen(nsrc)] 
  nfr = n_elements(frames)
  frames=reform(rebin(reform(frames,1,nfr),2,nfr),nfr*2) 
  ;tsize=[replicate(-1,n_elements(frames)), replicate(caldiam,n_elements(frames))]  
  tsize=[replicate(-1,n_elements(frames)/2), replicate(caldiam,n_elements(frames)/2)]  ;Fixed (?)
endif else begin
  frames = [src_start+indgen(nsrc)]
  nfr = n_elements(frames)
  frames=reform(rebin(reform(frames,1,nfr),2,nfr),nfr*2) 
  tsize=[replicate(-1,n_elements(frames))]
endelse
updown = reform(rebin([1,2],2,nfr),2*nfr)


; Modify variables to suit VAMPIRES LCVR switching
; (eventually this could be a switch).
; Different LCVR states are different cubes
frames = reform(rebin(reform(frames,1,nfr*2),2,nfr*2),nfr*4)
tsize = reform(rebin(reform(tsize,1,nfr*2),2,nfr*2),nfr*4)
updown = reform(rebin(reform(updown,1,nfr*2),2,nfr*2),nfr*4)
lcvrstate = reform(rebin([1,2],2,nfr*2),4*nfr)


;ddir_sky=ddir
set_plot,'x'
n_cubes=n_elements(frames)


; #### Currently copied from conica - need to tweak for VAMPIRES (eg
;                                     support LCVR)
; Construct cal4src matrix according to the cal4src flag
if(cal4src[0] eq -1) then begin              
   cal4src=intarr(n_cubes,n_cubes)      
   targ=where(tsize lt 0, complement=calib)
   if(calib[0] ne -1) then $
     for i=0,n_elements(targ)-1 do cal4src[calib,targ[i]]=1 
endif 
if(cal4src[0] eq -2) then begin
   cal4src = intarr(n_cubes,n_cubes)
   ngroups=nsrc/2
   for i=0,ngroups-1 do begin
     for j=0,3 do begin
       for k=j,nsrc*2-1,4 do begin
         cal4src[2*nsrc+k,4*i+j]=1
       endfor
     endfor
   endfor
endif
if(cal4src[0] eq -3) then begin
   cal4src = intarr(n_cubes,n_cubes)
   for i=0,2*nsrc-1 do begin
     cal4src[2*nsrc+i,i]=1
   endfor
endif
if(cal4src[0] eq -4) then begin
   cal4src = intarr(n_cubes,n_cubes)
   if cal_start GE 0 then begin
      for i=0,4*nsrc-1 do begin
         cal4src[4*nsrc+i,i]=1
      endfor
   endif
endif

; Define data structure to hold all the frame statistics:
ina=intarr(n_cubes)		; this is for the index of good frames (-1 implies bad frame)
f=fltarr(n_cubes,2)		; this holds frame statistic and err
fl= fltarr(n_cubes,10,2)          ; average and errors.  src#, aperture size, 0=mean/1=error
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}   ; ,good_frames:ina} 

; Get data structure to hold Header info and comments
skies = -1 ;No skies with VAMPIRES, but included for legacy code support
olog=make_olog_vamp(n_cubes,tsize,frames,skies)

; Begin to populate olog structure 
process_history=adate
olog.adate=adate
olog.rawdir=ddir
olog.comments=comments
olog.proc_hist=process_history
olog.cal4src=cal4src

if ( (size(discard_sigma))(0) eq 0) then begin
  discard_sigma=replicate(discard_sigma,5)
endif

supersky=fltarr(chipsize,chipsize) ; Don't do actual skies for emccd camera
flat=replicate(1.,chipsize,chipsize)

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s)
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cube_ix=0,n_cubes-1 do begin

    if tsize[cube_ix] LT 0 then begin
        prefix=src_prefix
    endif else begin
        prefix=cal_prefix
    endelse

      raw2cube_vampires,frames[cube_ix],supersky,datadir=ddir, 		$ ; Input params
	gain=flat,bad_pixels=bad_pixels,skies=skies,tsz=tsize[cube_ix], $
	cube=scube,headinfo=headinfo,fstats=fstats,flux=flux,   	$ ; Output Params
	noskyfit=noskyfit,plotall=plotall,setsquare=setsquare,          $
        prefix=prefix,extn = extn, dcube=dkcube,               $
        updown=updown[cube_ix],horiz_destripe=horiz_destripe,           $
        discard_last=discard_last, lcvrstate=lcvrstate[cube_ix],        $
        nocosmic=nocosmic, skyrad=skyrad, speckpos=speckpos, keepzeronum=keepzeronum

      info=size(scube)
      if(info[0] gt 2) then nspeck=info(3) else nspeck=1

      
      if (updown[cube_ix] lt 2) AND (lcvrstate[cube_ix] le 1) then $
        good_frames=flagbad_conica(fstats,scube,discard_sigma,tot_bad)
      good_index=where(good_frames ge 0)
      bad_frames=where(good_frames lt 0)


      ; %%%%% Plot out frame statistics
      if(identifier ne 'nosave') then $
        identifier = adate+"_"+string(frames[cube_ix],format="(I4.4)")+"a.ps"
      fnumstr=[string(frames[cube_ix],format="(I4.4)")]
      plota_nirc2,fstats,adate,headinfo,fnumstr,                        $
        process_history,comments,                                       $
        saturation_flag=saturation_flag,bad_frames=bad_frames,          $
        hardcopy=hardcopy,identifier=identifier

      ; Remove bad frame data from all arrays
      if ((size(bad_frames))[0] eq 0) then num_bad=0 else num_bad=n_elements(bad_frames)
      if(num_bad gt 0) then begin
          print,"Removing ",num_bad," flagged frames of bad data from this cube"
          scube=scube(*,*,good_index)
          fstats=fstats(*,*,good_index)
      endif


; %%%%% Reduce Aperture Photometry to averages 
; %%%%% and Populate framestats data structure.
  for flx=0,9 do framestats.phot[cube_ix,flx,0] =mean(flux[flx,*])
  for flx=0,9 do framestats.phot[cube_ix,flx,1] =stdev(flux[flx,*])
  framestats.xpk[cube_ix,*]        =[mean(fstats[0,*]),stdev(fstats[0,*])] 
  framestats.ypk[cube_ix,*]        =[mean(fstats[1,*]),stdev(fstats[1,*])] 
  framestats.totflx[cube_ix,*]     =[mean(fstats[2,*]),stdev(fstats[2,*])] 
  framestats.pkflx[cube_ix,*]      =[mean(fstats[3,*]),stdev(fstats[3,*])] 
  framestats.skybgr[cube_ix,*]     =[mean(fstats[4,*]),stdev(fstats[4,*])] 

; %%%%% Flux data to pass to diffcal, to go in csv file
  common datavals
  av_totalflux=mean(framestats.totflx[*,0])
  sd_totalflux=mean(framestats.totflx[*,1])
  av_pkflux=mean(framestats.pkflx[*,0])
  sd_pkflux=mean(framestats.pkflx[*,1])
  fluxvec=[av_totalflux,sd_totalflux,av_pkflux,sd_pkflux]

; %%%%% Populate olog data structure with defaults
  olog.instrument[cube_ix] =headinfo.instrument
  olog.nax1[cube_ix]       =headinfo.nax1
  olog.nax2[cube_ix]       =headinfo.nax2
  olog.t_int[cube_ix]      =headinfo.t_int
  olog.coadd[cube_ix]      =headinfo.coadd
  olog.filter[cube_ix]     =headinfo.filter
  olog.slit[cube_ix]       =headinfo.slit
  olog.optic_cfg[cube_ix]  =headinfo.optic_cfg
  olog.lyot[cube_ix]       =headinfo.lyot
  olog.grism[cube_ix]      =headinfo.grism
  olog.source_name[cube_ix]=headinfo.source_name
  olog.utc[cube_ix]        =headinfo.utc
  olog.date[cube_ix]       =headinfo.date
  olog.jd[cube_ix]         =headinfo.jd
  olog.elevation[cube_ix]  =headinfo.elevation
  olog.del_elev[cube_ix]   =headinfo.del_elev
  olog.airmass[cube_ix]    =headinfo.airmass
  olog.pa[cube_ix]         =headinfo.pa
  olog.del_pa[cube_ix]     =headinfo.del_pa
  olog.ra[cube_ix]         =headinfo.ra
  olog.dec[cube_ix]        =headinfo.dec
  olog.equinox[cube_ix]    =headinfo.equinox
  olog.mask                =headinfo.mask
  olog.flat_file           ='';flat_file

; Extra VAMPIRES-specific values:
  olog.emgain[cube_ix]     =headinfo.emgain
  olog.hwp[cube_ix]        =headinfo.hwp
  olog.timingpattern[cube_ix]=headinfo.timingpattern
  olog.imgRotAng[cube_ix]  =headinfo.imgRotAng
  olog.imgRotPad[cube_ix]  =headinfo.imgRotPad
  olog.imgRotPap[cube_ix]  =headinfo.imgRotPap
  olog.ADCstagepos[cube_ix]=headinfo.ADCstagepos
  olog.ADCp1angle[cube_ix] =headinfo.ADCp1angle
  olog.ADCp2angle[cube_ix] =headinfo.ADCp2angle
  olog.azimuth[cube_ix]    =headinfo.azimuth
  olog.localtime[cube_ix]  =headinfo.localtime












; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s) END MODULE 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


; %%%%% Now write out the cleaned source files as FITS
; NOTE: there is only the minimal fits header, and the file output does
;       not have a unique identifier: it will be overwritten by 
;       subsequent runs. This is deliberate - if you want to regenerate
;       a particular cube, run the script over again.
  filestring=[string(frames[cube_ix],format="(I4.4)")]
filestring='_'+prefix+filestring
if (updown[cube_ix] gt 0) then filestring=filestring+'_'+strtrim(updown[cube_ix],2)

if (lcvrstate[cube_ix] gt 0) then begin
    ; File will have suffixes '_1' or '_2' for up or down respectively,
    ; AND '_A' or '_B' for LCVR state 1 or state 2 respectively.
    if lcvrstate[cube_ix] eq 1 then lcvr_suffix='_A'
    if lcvrstate[cube_ix] eq 2 then lcvr_suffix='_B'
    filestring=filestring+lcvr_suffix
endif

  filename="cube"+filestring+".fits"
  olog.cube_fname[cube_ix,0]=filename
  olog.cube_fname[cube_ix,1]=filestring
  writefits,filename,scube
  if (n_elements(dcube) gt 1) then begin
    filename = "dcube"+filestring+".fits"
    olog.dk_fname[cube_ix] =  filename
    writefits,filename,dcube
  endif
  olog.cube_tsize[cube_ix]=tsize[cube_ix]
  olog.cube_sz[cube_ix,*]= (size(scube))[1:3]
endfor 
; END big loop over all cubes

; %%%%% Save the olog and framestats data structures for later use
save,olog,framestats,file='cubeinfo'+adate+'.idlvar'


end
