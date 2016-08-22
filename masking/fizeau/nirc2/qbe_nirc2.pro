; This was a Template that outputs cleaned NIRC2 data for further
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
;  skies    - File numbers of the sky files. Use -1 for dithered data.
; Optional inputs... but *either* fstart *or*
; frames and tsize must be set.
;  frames   - A vector of all used file numbers.
;  tsize    - A vector of the same length as frames, which has negative
;             sizes for targets and positive sizes for cals.
;  fstart   - A vector of starting file numbers, e.g. [732,743,754] 
;  nfs      - Number of data files corresponding to each element of
;             fstart. Use negative numbers for target stars, positive
;             numbers for calibrators e.g. [10,-8,10]. If not set,
;             then the blocks will be chopped up automatically,
;             assuming that there is no missing data.
;
; Optional Inputs:
;  discard_sigma - See below.

; Outputs: Files stored in IDL's current directory, with "cube" as 
;  part of the title. Make sure that the idl "cd" command is used 
;  before running this.
;
; Variant of cube_nirc2.pro which is used for new qball input 
;___________________________________________________________
;  Prime script information
;____________________________________________________________

pro qbe_nirc2, ddir, adate, flatfile, skies, frames=frames, tsize=tsize, fstart=fstart, nfs=nfs,$
    cal4src=cal4src, discard_sigma=discard_sigma, extn=extn, ddir_sky=ddir_sky,comments=comments,root_dir=root_dir

if not keyword_set(cal4src) then cal4src=-1
if not keyword_set(comments) then comments=''
if not keyword_set(extn) then extn =  '.fits.gz'
if keyword_set(fstart) then begin
    if not keyword_set(nfs) then begin
        nfs = fstart[1:*]-fstart
        if min(nfs) le 0 then begin
            print, 'ERROR: fstart must increase monotonically if you don''t set nfs'
            stop
        endif
        endfile = max(fstart)+1
        ;;maxendfile is the last file that we'll consider to be part
        ;;of the set.
        maxendfile = max(fstart) + nfs(n_elements(nfs)-1)
        while (file_test(ddir +'/n' + string(endfile, format='(I04)') + extn) and endfile lt maxendfile) $
          do endfile++
        nfs = -[nfs,endfile-max(fstart)]
    endif else if (n_elements(fstart) ne n_elements(nfs)) then begin
        print, 'ERROR: fstart and nfs must have the same number of frames'
        stop
    endif 
    frames = make_frame_vect(fstart,  nfs,  tsize = tsize)
endif else if not (keyword_set(frames) and keyword_set(tsize))  then begin
   print, 'ERROR: You must either set frames and tsize or set fstart and nfs'
   stop
endif else if (n_elements(frames) ne n_elements(tsize)) then begin
    print, 'ERROR: frames and tsize must have the same number of frames'
    stop
endif

; Restore Flat.
restore, flatfile

if not keyword_set(ddir_sky) then ddir_sky=ddir
;;Make sure that ddir and ddir_sky finish with a '/' !!!TODO


;____________________________________________________________
;  Analysis Options
;____________________________________________________________
noskyfit=0	       ; 0 = Chip-periphery sky subtraction [FOR FAINT SOURCES!]
		       ; 1 = NO Sky Subtraction (Recommended for Golay/subframe)
                       ; 2 = Don't even compute sky bgr! (=1 still puts it in stats)
setsquare=1	       ; 0 = Preserve (non-square) aspect
		       ; 1 = Trim to square (auto)
                       ; 256 = trim to array size 256 pix
if not keyword_set(discard_sigma) then $
  discard_sigma=[-2., -2., 3.5,  3.5, 0.4]
                       ; 0 = manual discard (examine frames on-screen)
                       ; -1 = No discard procedure
                       ; Vector [3.,3.,3.0,5.0] Removes Speckle frames whose statistics 
		       ;  are greaters than (discard_sigma) standard deviations from the mean
                       ;  Any negative elements means ignore testing that variable (but take
                       ;  care with setting the first element to -1 or *ALL* testing is abandoned)
		       ; **NOTE** if discard_sigma is a 5 element 
		       ; vector, then apply different sigma cuts for
		       ; each diagnostic variable.  
		       ; USAGE: discard_sigma=
                       ;   [x_pos_sigma,y_pos_sigma,total_counts_sigma,
                       ;       peak_pixel_sigma, fraction of peak median]
                       ;discard_sigma[4] is designed to reject low-strehl frames
                       ;For dithered data, you can set e.g. x_pos_sigma to 1.5.
save_dcube = 1         ;Save the a cube of darks (for power spectrum analysis)
;_____________________________________________________________
;  Information for Plotting
;_____________________________________________________________
plotall =1             ; 0 = no plots to screen
                       ; 1 = plot the images only. Not implemented...
                       ; 2 = Plot the sky subtraction and the
                       ;  images on top. 
hardcopy=0	       ; 0 = no hardcopy output
		       ; 1 = print all output
identifier='save'      ; save ps plots in data directory (default)
;identifier='nosave'      ; save ps plots in data directory (default)
                       ; 'nosave'= do not save output plots
;_____________________________________________________________
; The rest is done automatically.
;_____________________________________________________________
process_history=adate
if (noskyfit eq 0) then process_history=process_history+"Sky Subtraction;" else $
                     process_history=process_history+"No Sky Subtraction;"

set_plot,'x' ; bug fix - prevents runaway printing in cases after a crash
             ;   (e.g. all 100 sky subtract histograms can be printed out)

nframes=n_elements(frames)
scsc=uniq(tsize)           ; This tells how many cubes need to be made
n_cubes=n_elements(scsc)

if(cal4src[0] eq -1) then begin                ; Default is to assume all cals will be used
   cal4src =  intarr(n_cubes,  n_cubes)        ; for each source. Set this var in flagging later.
   targ  =  where(tsize(scsc) lt 0,  complement = calib)
   if(calib[0] ne -1) then $
     for i =  0, n_elements(targ)-1 do cal4src[calib,targ[i]] =  1  
endif
scsc=[-1,scsc]              ; File number intervals start at zero (helps later)

; Define data structure to hold all the frame statistics:
ina=intarr(n_cubes)		; this is for the index of good frames (-1 implies bad frame)
f=fltarr(n_cubes,2)		; this holds frame statistic and err
fl= fltarr(n_cubes,10,2)          ; average and errors.  src#, aperture size, 0=mean/1=error
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}   ; ,good_frames:ina} 


; Get data structure to hold Header info and comments
olog=make_olog(n_cubes,tsize,frames,skies)

; Begin to populate olog structure 
olog.adate=adate
olog.rawdir=ddir
olog.comments=comments
olog.proc_hist=process_history
olog.cal4src=cal4src
olog.flat_file=flatfile

if ( (size(discard_sigma))(0) eq 0) then begin
  discard_sigma=replicate(discard_sigma,5)
endif

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s)
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cube_ix=0, n_cubes-1 do begin
  print,'Up to cube_ix = ',cube_ix
  scframes=frames[scsc[cube_ix]+1:scsc[cube_ix+1]]    ; filenums in this cube

  ;;Setsquare: MJI isn't sure what this really does. 
  if(setsquare eq 0) then print,'ERROR ### Cannot preserve full frame for dithered data'
  if(setsquare eq 1) then setsquare=256
  raw2cube_dither_nirc2,scframes,datadir=ddir, 		$ ; Input params
	gain=flat,bad_pixels=bad_pixels,                        $
	cube=scube,headinfo=headinfo,fstats=fstats,flux=flux,	$ ; Output Params
	noskyfit=noskyfit,setsquare=setsquare,plotall=plotall,  $
        dcube = dcube,  skies = skies, extn=extn, skydir=ddir_sky, /medsub_ssky,root_dir=root_dir;;,  /destripe
  if (save_dcube ne 1) then dcube = -1

  info=size(scube)
  if(info[0] gt 2) then nspeck=info(3) else nspeck=1

good_frames=flagbad_nirc2(scube,fstats,discard_sigma,scframes)
good_index=where(good_frames ge 0)
bad_frames=where(good_frames lt 0)

; %%%%% Plot out frame statistics
  if(identifier ne 'nosave') then $
  identifier = strcompress(adate+"_"+string(scframes(0))+"-"+ $
                 string(scframes(nspeck-1))+"a",/remove_all)+".ps"
  nsrc=n_elements(scframes)
  fnumstr=[string(scframes[0],format="(I4.4)")+'-'+string(scframes[nsrc-1],format="(I4.4)")]
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
  olog.mask=inquire('mask',  olog)


; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Target Star(s) END MODULE 
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print,''

; %%%%% Now write out the cleaned source files as FITS
; NOTE: there is only the minimal fits header, and the file output does
;       not have a unique identifier: it will be overwritten by 
;       subsequent runs. This is deliberate - if you want to regenerate
;       a particular cube, run the script over again.
  nsrc=n_elements(scframes)
  filestring=[string(scframes[0],format="(I4.4)")]
  filename="cube"+filestring+".fits"
  olog.cube_fname[cube_ix,0]=filename
  olog.cube_fname[cube_ix,1]=filestring
  writefits,filename,scube
  if (n_elements(dcube) gt 1) then begin
    filename = "dcube"+filestring+".fits"
    olog.dk_fname[cube_ix] =  filename
    writefits,filename,dcube
  endif
  olog.cube_tsize[cube_ix]=tsize[scsc[cube_ix+1]]
  olog.cube_sz[cube_ix,*]= (size(scube))[1:3]
endfor 
; END big loop over all cubes

; %%%%% Save the olog and framestats data structures for later use
save,olog,framestats,file='cubeinfo'+adate+'.idlvar'


end

