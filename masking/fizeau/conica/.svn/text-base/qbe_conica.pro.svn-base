;######## RELEASE CLIP HERE ######## 
; This file is built from a merger between qbe_nirc2 and cube_conica
; This was a Template that outputs cleaned CONICA data for further
; analysis later. Now it is a procedure. 
; Src and Cal input files are to be glued into tw
; Philosophy is to do all basic image-plane processing (eg sky
; subtract etc) but leave the FFT part for later routines.
; Also generate diagnostics, and reject abberrant frames.
; 
; Inputs:
;  ddir     - Full path of data directory
;  adate    - Analysis date or label for output files.
;  flatfile - idl variable file containing a flat map and bad pixels.
;  fstart   - A vector of starting file numbers, e.g. [732,743,754] 
;  nfs      - Number of data files corresponding to each element of
;             fstart. Use negative numbers for target stars, positive
;             numbers for calibrators e.g. [10,-8,10]
;  skies    - File numbers of the sky files. Use -1 for dithered data.
;
;
; Optional Inputs:
;  discard_sigma - See below.
;  update_flat - Gives the option to update the bad pixel map in the
;                flat file with the bad pixels found in getskyconica
;  discard_ends - Ignores the first and last frame of each cube, since
;                 the last is a sum of the rest of the cube, and the
;                 first suffers from reset errors (from NACO manual).
;                 THIS IS SET BY DEFAULT. To discard only the last frame,
;                 set discard_ends=2. To keep all frames, set discard_ends=3
;  noplots      - Skips all plotting for all subprograms, so that the
;                 pipeline can be run as a batch process. (same as plotall?
;
; Outputs: Files stored in IDL's current directory, with "cube" as 
;  part of the title. Make sure that the idl "cd" command is used 
;  before running this.
;
; Variant of cube_conica.pro which is used for new qball input 
;___________________________________________________________
;  Prime script information
;____________________________________________________________


pro qbe_conica, ddir, adate, flatfile, frames, skies, tsize,   	 	$
                discard_sigma = discard_sigma, prefix = prefix, extn = extn, 	   			$
	        			updown = updown, horiz_destripe = horiz_destripe, noskyfit = noskyfit,$
	        			discard_ends = discard_ends, ddir_sky = ddir_sky, plotall = plotall,	$
	    					setsquare = setsquare, medsub = medsub, hardcopy = hardcopy,  		  	$
								ps_save = ps_save, pause=pause,update_flat=update_flat,noplots=noplots

restore, flatfile; Restore Flat.
;frames = make_frame_vect(fstart,  nfs,  tsize = tsize)
; 
comments='OK'

if not keyword_set(ddir_sky) or n_elements(ddir_sky) eq 0 then ddir_sky=ddir
if not keyword_set(medsub) then medsub=0
if not keyword_set(pause) then pause=0.1
if not keyword_set(discard_ends) then discard_ends=1
;____________________________________________________________
;  Analysis Options
;____________________________________________________________
if not keyword_set(noskyfit) then noskyfit=0
    	    	        ; 0 = Chip-periphery sky subtraction [FOR FAINT SOURCES!]
		        				; 1 = NO Sky Subtraction (Recommended for Golay/subframe)
                    ; 2 = Don't even compute sky bgr! (=1 still puts it in stats)

if not keyword_set(setsquare) then setsquare=1
    	    	        ; 0 = Preserve (non-square) aspect
		        				; 1 = Trim to square (auto)
                    ; 256 = trim to array size 256 pix
if not keyword_set(discard_sigma) then discard_sigma=[-2., -2., 3.5,  3.5, 0.4]
                    ; 0 = manual discard (examine frames on-screen)
                    ; -1 = No discard procedure
                    ; Vector [3.,3.,3.0,5.0] Removes Speckle frames whose statistics 
						        ;  are greaters than (discard_sigma) standard deviations from the mean
                    ;  Any negative elements means ignore testing that variable (but take
                    ;  care with setting the first element to -1 or *ALL* testing is abandoned)
		        				; **NOTE** if discard_sigma is a 4 element 
		        				; vector, then apply different sigma cuts for
		        				; each diagnostic variable.  
		        				; USAGE: discard_sigma=
                    ;   [x_pos_sigma,y_pos_sigma,total_counts_sigma,
                    ;       peak_pixel_sigma, fraction of peak median]
                    ;discard_sigma[4] is designed to reject low-strehl frames
                    ;For dithered data, you can set e.g. x_pos_sigma to 1.5.
save_dcube = 1      ;Save the a cube of darks (for power spectrum analysis)
;_____________________________________________________________
;  Information for Plotting
;_____________________________________________________________
if not keyword_set(plotall) then plotall = 0
    	               ; 1 = no plots to screen
                     ; 0 = plot to screen
plotall = ~plotall   ; invert so it works further down

if not keyword_set(hardcopy) then hardcopy=0
									; 0 = no hardcopy output
		        			; 1 = print all output
if not keyword_set(ps_save) then ps_save = 0
if ps_save eq 1 then identifier='save' else identifier='nosave'
  								; save ps plots in data directory (default)
                  ; 'nosave'= do not save output plots
;_____________________________________________________________
; The rest is done automatically.
;_____________________________________________________________
process_history=adate
if (noskyfit eq 0)  then process_history=process_history+"Sky Subtraction;" $
										else process_history=process_history+"No Sky Subtraction;"
if(setsquare eq 0) then print,'ERROR ### Cannot preserve full frame for dithered data'
;setquare=1 => auto, sets setsquare to half the size of the array
if(setsquare eq 1) then setsquare=fix((size(flat))(1)/2)
skymask=fix(setsquare/2.0)

; Sort out file prefixes ... if we only have one passed, then replicate it...
if (n_elements(prefix) le 1) then prefix = replicate(prefix, n_elements(frames))
	
    ;build skies
    if(n_elements(skies) gt 1) then begin 
    	skyextn = skies[0]
	skies = skies[1:n_elements(skies)-1]
	skyprefix = skies[n_elements(skies)/2:n_elements(skies)-1]
	skies = skies[0:(n_elements(skies)/2)-1]
	nskies=n_elements(skies) 
	for i=0,nskies-1 do begin	; BEGIN LOOP OVER SKY FILES
	    cleansky_conica,skies[i],prefix=skyprefix[i],extn=skyextn,bad_pixels=bad_pixels,  $
            datadir=ddir_sky,root_dir=root_dir,supersky=ssky
     	    if(i eq 0) then skycube=ssky else skycube=[[[skycube]],[[ssky]]]
    	endfor
    	if(nskies gt 1) then begin
     	    supersky=total(skycube,3)/nskies
     	    for i=0,nskies-1 do begin
            	skycube[*,*,i]=skycube[*,*,i]-supersky
      	    	skycube[*,*,i]=skycube[*,*,i]/flat
     	    endfor
    	endif else begin
    	    supersky=skycube
    	    skycube=skycube/flat
    	endelse
    endif else begin  ; OK so we have to make up a sky with all the images we have ...
    	if(skies eq "-1") then begin   ; case: make up a single sky frame for all data.
    	    skies = fix(skies)
	    nskies=1
    	    getssky_conica,frames,supersky,datadir=ddir,prefix=prefix,extn = extn, gain=flat, $
              starmask=skymask,handsky=dohandsky, sky_lr=sky_lr, bad_pixels=bad_pixels, 				$
              polzflag=polzflag, medsub=medsub,new_bad_pixels=new_bad_pixels, $
              discard_ends=discard_ends
        endif 
    	if(skies eq "-2") then begin   ; case: make up 2 sky frames for src and cal
    	    skies = fix(skies)
	    nskies=2
    	    fsrc=where(tsize lt 0) 
            getssky_conica,frames(fsrc),supersky_src,datadir=ddir,prefix=prefix[fsrc],extn = extn,gain=flat, $
              starmask=skymask,handsky=dohandsky, sky_lr=sky_lr, bad_pixels=bad_pixels,polzflag=polzflag,$
              medsub=medsub,new_bad_pixels=new_bad_pixels_src,discard_ends=discard_ends
    	    fcal=where(tsize ge 0)
            getssky_conica,frames(fcal),supersky_cal,datadir=ddir,prefix=prefix[fcal],extn = extn,gain=flat, $ 
              starmask=skymask,handsky=dohandsky, sky_lr=sky_lr, bad_pixels=bad_pixels,polzflag=polzflag,$
              medsub=medsub,new_bad_pixels=new_bad_pixels_cal,discard_ends=discard_ends
            new_bad_pixels=new_bad_pixels_src+new_bad_pixels_cal
            supersky=[[[supersky_src]],[[supersky_cal]]]
        endif
        if (skies eq "-3") then begin ;;separate dithers for each file (WORK IN PROGRESS)
            blocks=[0,uniq(tsize)+1];;this will make one more entry in blocks than there are actual blocks
            skies=fix(skies)
            nskies=n_elements(blocks)-1
            for skynum=0,nskies-1 do begin
                getssky_conica,frames(blocks[skynum]:blocks[skynum+1]-1),supersky_temp,datadir=ddir,prefix=prefix[blocks[skynum]:blocks[skynum+1]-1],extn = extn,gain=flat, $
                  starmask=skymask,handsky=dohandsky, sky_lr=sky_lr, bad_pixels=bad_pixels,polzflag=polzflag,$
                  medsub=medsub,new_bad_pixels=new_bad_pixels,discard_ends=discard_ends
                if skynum gt 0 then supersky=[[[supersky]],[[supersky_temp]]] else supersky=supersky_temp
            endfor
        endif
    endelse
;Write out the ssky to view later...
    save,supersky,file='supersky_saved.idlvar'
   
;;update the bad pixel map from the flat:
if keyword_set(update_flat) then begin
    bp_map=0*new_bad_pixels
    bp_map[bad_pixels]=1
    all_bad=bp_map+new_bad_pixels;;old bad pixel map
    bad_pixels=where(all_bad gt 0.5) ;;worse more than 50% of the time
    print,'save updated bad pixel map? Use  save,flat,bad_pixels,filename=flatfile  .c to continue'
    stop ;; check that it worked, .c to continue
    ;;save,flat,bad_pixels,filename=flatfile
endif

set_plot,'x' ; bug fix - prevents runaway printing in cases after a crash
             ;   (e.g. all 100 sky subtract histograms can be printed out)

;nframes=n_elements(frames)
; scsc=uniq(tsize)           ; This tells how many cubes need to be made
;n_cubes=n_elements(scsc)
n_cubes=n_elements(frames)

if not keyword_set(updown) then updown=replicate(0,n_cubes)

cal4src =  intarr(n_cubes,  n_cubes) ; for each source. Set this var in flagging later.
targ  =  where(tsize lt 0,  complement = calib)
if(calib[0] ne -1 and targ[0] ne -1) then $
  for i =  0, n_elements(targ)-1 do cal4src[calib,targ[i]] =  1  

;scsc=[-1,scsc]              ; File number intervals start at zero (helps later)

; Define data structure to hold all the frame statistics:
ina=intarr(n_cubes)		; this is for the index of good frames (-1 implies bad frame)
f=fltarr(n_cubes,2)		; this holds frame statistic and err
fl= fltarr(n_cubes,10,2)          ; average and errors.  src#, aperture size, 0=mean/1=error
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}   ; ,good_frames:ina} 

; Get data structure to hold Header info and comments
olog=make_olog(n_cubes,tsize,frames,skies)

;Begin to populate olog structure 
olog.adate=adate
olog.rawdir=ddir
olog.comments=comments
olog.proc_hist=process_history
olog.cal4src=cal4src

if ( (size(discard_sigma))(0) eq 0) then begin
  discard_sigma=replicate(discard_sigma,5)
endif

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s)
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cube_ix=0, n_cubes-1 do begin
    ;;rather than change raw2cube_conica, it is easiest to add in an
    ;;extra few lines here to pass the correct supersky in the case of
    ;;separate skies for each observing block (skies[0] eq 3)
    if fix(skies[0]) eq -3 then begin
        this_block=(where(blocks gt cube_ix))[0]-1
        current_supersky=supersky[*,*,this_block]
    endif else current_supersky=supersky
  ;scframes=frames[scsc[cube_ix]+1:scsc[cube_ix+1]]    ; filenums in this cube
  if(setsquare eq 0) then print,'ERROR ### Cannot preserve full frame for dithered data'
  ;setquare=1 => auto, sets setsquare to half the size of the array
	if(setsquare eq 1) then setsquare=fix((size(flat))(1)/2)
        raw2cube_conica,frames[cube_ix],current_supersky,datadir=ddir, 	            $ ; Input params
 	    gain=flat,bad_pixels=bad_pixels, pause=pause,tsz=tsize(cube_ix),$
	    cube=scube,headinfo=headinfo,fstats=fstats,flux=flux,	$ ; Output Params
	    noskyfit=noskyfit,setsquare=setsquare,plotall=plotall,    	$
  	    dcube = dcube,  skies = skies, extn=extn, quad=quad,  	$
	    prefix=prefix[cube_ix], updown=updown[cube_ix],     	$
	    horiz_destripe=horiz_destripe, discard_ends=discard_ends
  if (save_dcube ne 1) then dcube = -1
  info=size(scube)
  if(info[0] gt 2) then nspeck=info(3) else nspeck=1

;;good_frames=flagbad_conica(scube,fstats,discard_sigma,scframes[cube_ix])
good_frames=flagbad_conica(fstats,scube,discard_sigma,tot_bad)
good_index=where(good_frames ge 0)
bad_frames=where(good_frames lt 0)
; %%%%% Plot out frame statistics
  if(identifier ne 'nosave') then $
;;  identifier = strcompress(adate+"_"+string(scframes(0))+"-"+ $
;;                 string(scframes(n_elements(scframes)-1))+"",/remove_all)+".ps"
;;  nsrc=n_elements(scframes)
;;  fnumstr=[string(scframes[0],format="(I4.4)")+'-'+string(scframes[nsrc-1],format="(I4.4)")]
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

; %%%%% Populate olog data structure.
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
  olog.mask                =headinfo.mask;inquire('mask',  olog)
  olog.quad[cube_ix]       =quad

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s) END MODULE 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print,''
; %%%%% Now write out the cleaned source files as FITS
; NOTE: there is only the minimal fits header, and the file output does
;       not have a unique identifier: it will be overwritten by 
;       subsequent runs. This is deliberate - if you want to regenerate
;       a particular cube, run the script over again.
;;  nsrc=n_elements(scframes)
;;  filestring=[string(scframes[0],format="(I4.4)")]

  filestring=[string(frames[cube_ix],format="(I4.4)")]
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


;print, olog.cal4src

; Generate a single postscript file from all "a.ps" files for ease of clockthru
; (optional this can be commented out)
spawn, 'psjoin *a.ps > aplots.ps'

; %%%%% Save the olog and framestats data structures for later use
save,olog,framestats,file='cubeinfo'+adate+'.idlvar'


end

