;$Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $
;$Log: calc_bispect.pro,v $
;Revision 1.25  2010/06/22 06:20:15  mireland
;Pretty sure this was very minor...
;
;Revision 1.24  2010/04/08 20:47:02  dbernat
;switched ph_all to cvis_all
;
;Revision 1.23  2010/03/26 11:44:29  mireland
;Added a !ROOT_DIR system variable, changed the inputs for qbe_nirc2.pro
;so that you can specify all indices or just starting indices, and fixed a bug
;where bad data flaging was turned off.
;
;Revision 1.22  2010/03/24 03:42:42  snert
;Added clip marker for release version
;
;Revision 1.21  2010/03/16 01:42:51  mireland
;Added a cp output, so you dont need to atan(bs, /phase)
;
;Revision 1.20  2010/03/16 01:30:40  mireland
;Added v2_sig and cp_sig to calc_bispect.
;
;Revision 1.19  2009/09/30 02:45:15  snert
;Moved the plot command outside the for loop so that the plot window doesn't come to the front of the screen every time a new power spectrum is plotted.
;
;Revision 1.18  2009/04/12 11:30:29  snert
;Not quite sure what these changes are...
;
;Revision 1.17  2009/04/01 00:04:42  mireland
;Added a coherent averaging keyword to calc_bispect
;
;Revision 1.16  2009/03/03 07:11:11  snert
;MJI commiting stuff for Conica.
;
;Revision 1.15  2008/11/02 21:50:35  mireland
;A couple of bug fixes inspired by Woody.
;
;Revision 1.14  2008/10/22 06:22:47  mireland
;calibrate_v2_cp had a conflict and isn't tested. But I think that there are
;basically bugfixes only here.
;
;Revision 1.13  2008/09/17 18:00:40  dbernat
;Baseline Phase information now saved (ph_all)
;
;Revision 1.12  2008/09/09 20:12:14  dbernat
;now saves all bs and v2
;
;Revision 1.11  2008/05/21 23:11:28  mireland
;Not quite sure what all these changes are - but now they are commmited anyway.
;
;Revision 1.10  2007/10/06 09:44:37  mireland
;Various changes, mostly relating to the re-analysis of USco data.
;
;Revision 1.9  2007/06/18 17:19:13  mireland
;Bugfixes...
;
;Revision 1.8  2007/06/15 00:31:28  mireland
;Still working on the covariance matrix stuff...
;
;Revision 1.7  2007/06/12 23:52:40  mireland
;Lots of changes: improved cosmic ray rejection, better naming and
;directory stuff, modifications to closure-phase histogram plots...
;
;Revision 1.6  2007/05/18 19:11:10  mireland
;Changed some defaults (histogram cutoff, window size) and removed a couple of annoying
;plots that flashed up.
;
;Revision 1.1  2006/01/12 23:50:23  mireland
;Commit of the new .pro calc_bispect and calibrate scripts, and the
;first LWS code version.
;
;Revision 1.3  2005/12/21 21:08:33  mireland
;Found a bug: there was no multiplication by the window function. Now fixed.
;
;Revision 1.2  2005/12/20 21:51:55  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and $Log: calc_bispect.pro,v $
;Added $Id: calc_bispect.pro,v 1.24 2010/04/08 20:47:02 dbernat Exp $ and Revision 1.25  2010/06/22 06:20:15  mireland
;Added $Id: calc_bispect.pro,v 1.24 2010/04/08 20:47:02 dbernat Exp $ and Pretty sure this was very minor...
;Added $Id: calc_bispect.pro,v 1.24 2010/04/08 20:47:02 dbernat Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.24  2010/04/08 20:47:02  dbernat
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and switched ph_all to cvis_all
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.23  2010/03/26 11:44:29  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Added a !ROOT_DIR system variable, changed the inputs for qbe_nirc2.pro
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and so that you can specify all indices or just starting indices, and fixed a bug
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and where bad data flagging was turned off.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.22  2010/03/24 03:42:42  snert
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Added clip marker for release version
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.21  2010/03/16 01:42:51  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Added a cp output, so you dont need to atan(bs, /phase)
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.20  2010/03/16 01:30:40  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Added v2_sig and cp_sig to calc_bispect.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.19  2009/09/30 02:45:15  snert
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Moved the plot command outside the for loop so that the plot window doesn't come to the front of the screen every time a new power spectrum is plotted.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.18  2009/04/12 11:30:29  snert
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Not quite sure what these changes are...
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.17  2009/04/01 00:04:42  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Added a coherent averaging keyword to calc_bispect
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.16  2009/03/03 07:11:11  snert
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and MJI commiting stuff for Conica.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.15  2008/11/02 21:50:35  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and A couple of bug fixes inspired by Woody.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.14  2008/10/22 06:22:47  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and calibrate_v2_cp had a conflict and isn't tested. But I think that there are
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and basically bugfixes only here.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.13  2008/09/17 18:00:40  dbernat
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Baseline Phase information now saved (ph_all)
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.12  2008/09/09 20:12:14  dbernat
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and now saves all bs and v2
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.11  2008/05/21 23:11:28  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Not quite sure what all these changes are - but now they are commmited anyway.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.10  2007/10/06 09:44:37  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Various changes, mostly relating to the re-analysis of USco data.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.9  2007/06/18 17:19:13  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Bugfixes...
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.8  2007/06/15 00:31:28  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Still working on the covariance matrix stuff...
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.7  2007/06/12 23:52:40  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Lots of changes: improved cosmic ray rejection, better naming and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and directory stuff, modifications to closure-phase histogram plots...
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.6  2007/05/18 19:11:10  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Changed some defaults (histogram cutoff, window size) and removed a couple of annoying
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and plots that flashed up.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.1  2006/01/12 23:50:23  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Commit of the new .pro calc_bispect and calibrate scripts, and the
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and first LWS code version.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Revision 1.3  2005/12/21 21:08:33  mireland
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and Found a bug: there was no multiplication by the window function. Now fixed.
;Added $Id: calc_bispect.pro,v 1.25 2010/06/22 06:20:15 mireland Exp $ and to important files, and added the g18 mask
;for PHARO's inquire.
;
;######## RELEASE CLIP HERE ######## 
;This procedure takes a bunch of data cubes in units of
;photons and outputs idl variable files containing their 
;bispectra and covariance matrices.
;
;NB Stuff this procedure doesn't do yet:
; Average splodges in bispectrum/power space
; Variables:
; cubeinfo_file The cubeinfo file, that contains an olog. 
; Options:
;  root_dir:    The root directory for templates (e.g. ~/code for
;               Mike's directory structure or ~/code/masking for 
;               Tom's directory structure). If this isn't set, you
;               need to devsysv, '!ROOT_DIR', 'THE_DIRECTORY'
;  ddir:        Set this to the location of the cubes, if it isn't
;               the same as the cubeinfo file.
; reset:        Set this keyword if the data had already been
;               processed with some special options, but you want to
;               re-process with default options.
; plotall 	       ; 0 = no hardcopy output (default)
;		       ; 1 = plot all output
; 	               ; 2 = plot final result but nothing else
; nosave        Set this keyword to
; coherent_vis  Do we coherently integrate vis over a splodge??
;               (default 1=yes)
;  coherent_bs  Do we coherently integrate the bispectrum over a
;               splodge??? (default 1=yes)
; window_type   1 for Super-Gaussian, 0 for none  (default 1)
; window_size   See below - the size of the window is generally set by
;               the hole size and wavelength.
; median_sub    For data with no darks or strong sky bg, will
;               subtract...   
;                =1 median of the data outside of the window.
;                =2 mean of the data outside of the window
;                (non-gaussian bg).
; disp_windowing: displays both windowed and unwindowed data, useful for
;                 checking that you are windowing out companions etc

pro calc_bispect, cubeinfo_file,  root_dir=root_dir, ddir = ddir,  reset = reset,  plotall =plotall,  nosave = nosave, $
                  coherent_vis = coherent_vis,  coherent_bs = coherent_bs,  window_type = window_type, add_noise = add_noise, $
                  window_size = window_size,  mf_file = mf_file,  median_sub = median_sub,  n_blocks = n_blocks, $
                  subtract_bs_bias = subtract_bs_bias, window_mult=window_mult, special=special, cav=cav, rebin=rebin,disp_windowing=disp_windowing,all_ps=all_ps

; File and Directory Options
showps = 0
;; Now set the default root_dir if we can
defsysv, '!ROOT_DIR', exists=exists
if exists then root_dir=!ROOT_DIR
if not keyword_set(root_dir) then begin
 print, 'ERROR: Must set root_dir keyword or !ROOT_DIR system variable.'
endif
;;"special" is a special string to tack on to the bs file
if not keyword_set(special) then special=''
if not keyword_set(cav) then cav=1

;The analysis dir should be where the cubeinfo file lives.
;The cube directory can be specified in olog: if not, it is assumed to
;be the same as the cubeinfo file's directory.
analysis_dir =  ''
if (keyword_set(ddir) eq 0) then ddir =  ''
pos = strpos(cubeinfo_file,'/',  /reverse_search)
if (pos ne -1) then begin
    analysis_dir =  strmid(cubeinfo_file,  0,  pos) + '/'
    if (keyword_set(ddir) eq 0) then ddir =  analysis_dir
    cubeinfo_file = strmid(cubeinfo_file,  pos+1)
endif
restore, analysis_dir +cubeinfo_file 
if not keyword_set(ddir) then $
  if (olog.cubedir ne '') then ddir =  olog.cubedir

if (keyword_set(reset)) then olog.logflag = 0
;if (keyword_set(ddir) eq 0) then plotall = 0
if (keyword_set(plotall) eq 0) then plotall = 0
if (keyword_set(median_sub) eq 0) then median_sub = 0
if (keyword_set(subtract_bs_bias) eq 0) then subtract_bs_bias = 0
                                ;These next options shouldn't be changed...
in_prefix = 'cube'
in_ext = '.fits'
print_prefix = 'test'
print_ext = '.ps'
;____________________________________________________________
;  Analysis Options
;____________________________________________________________

if (keyword_set(mf_file) eq 0) then begin 
    mf_file =  inquire('template',olog,ix=0)
    if (olog.logflag ne 0) then mf_file = plog.mf_file
    print,  'Using Matched Filter file: ',  mf_file

endif
if (keyword_set(n_blocks) eq 0) then begin
    n_blocks =  inquire('n_blocks',  olog,  ix = 0)
endif
closing_tri_pix=0  ; set this to zero as default for MJI-style mf files
restore,  root_dir + '/templates/' + mf_file
;restore, root_dir + mf_file

;Calculate a default image size
imsize =  filter[0]/hole_diam/rad_pixel
;imsize =  10 ;\lambda/hole_diam in pixels

if (keyword_set(window_size) eq 0) then begin
    window_size = 0.8*imsize ;1.3 suppresses, but still includes thr first airy ring.
                                ; 0.8 is slightly larger than S/N optimal.
                                ; This is the radius from the image centroid
                                ; that the window goes down to half it's peak
                                ; value.
    if (keyword_set(window_mult) eq 0) then begin
        print,  'window_size set to the default of: ',  window_size
    endif else begin
        window_size *= window_mult
        print,  'window_size set to ', window_size,' , a fixed multiple of the default.'
    endelse
endif else print,  'manual window_size change to: ',  window_size

if (keyword_set(coherent_vis) eq 0) then coherent_vis =  1 
if (keyword_set(coherent_bs) eq 0) then coherent_bs  =  1 
if (keyword_set(window_type) eq 0) then window_type =  1  
;_____________________________________________________________
; The rest is done automatically.
;

set_plot,'x' 

n_runs=(size(olog.cube_fname))[1]
v2 = fltarr(n_baselines)
v2_cov = fltarr(n_baselines,n_baselines)
bs = complexarr(n_bispect)
bs_var = fltarr(2,n_bispect)
cp_var = fltarr(n_bispect)
bs_cov = fltarr(2,n_cov)
bs_v2_cov = fltarr(n_baselines,n_holes-2)

; Begin analysis.

;if plotall eq 0 then window,  0,  xsize =1100,  ysize = 500 

loadct,  39
for i=0, n_runs-1 do begin

   if (plotall ne 0) then begin
       set_plot, 'ps'
       device,ysize=14,xsize=21,bits=8,/color,yoffset=0,xoffset=0;,/landscape
   endif ;else begin
       ;window,  0,  xsize =1100,  ysize = 500 
   ;endelse

    filename = olog.cube_fname[i,0]
    cube = readfits(ddir+filename)
    dimx = (size(cube))[1]
    dimy = (size(cube))[2]

    if i eq 0 then all_ps=fltarr(dimx,dimy,n_runs)

    if ((size(cube))[0] eq 2) then n_frames = 1 else n_frames = (size(cube))[3]
    if (keyword_set(rebin)) then begin
        n_frames = n_frames/2
        cube = cube[*,*,0:2*(n_frames-1):2] + cube[*,*,1:2*(n_frames-1)+1:2]
    endif
    if (keyword_set(cav)) then begin
        n_frames = n_frames/cav
        cube = cube[*,*,0:n_frames*cav-1]
        cube = rebin(cube, dimx, dimy, n_frames)
    endif
    if (keyword_set(add_noise)) then cube += add_noise*randomn(seed, dimx, dimy, n_frames)
    dim_min =  min([dimx, dimy])
    window =  fltarr(dimx, dimy)
    if (window_type eq 1) then $
      window[0:dim_min-1, 0:dim_min-1] = exp(-(dist(dim_min)/window_size*0.91244)^4) $ 
    else if (window_type eq 0) then $
      window[*] =  1.0
    if(median_sub ne 0) then begin
        if (window_size gt dimx/2) then stop
        if (dimx ne dimy) then stop ;Following line will not work
        w =  where(shift(dist(dimx), dimx/2, dimx/2) gt 1.2*window_size)
        if median_sub eq 1 then for k = 0, n_frames-1 do cube[*, *, k] -= median((cube[*, *, k])[w])
        if median_sub eq 2 then for k = 0, n_frames-1 do cube[*, *, k] -= mean((cube[*, *, k])[w]) 
    endif
    ft_cube = complexarr(dimx,dimy,n_frames)
    if (olog.dk_fname[i] ne '') then begin
        dcube = readfits(ddir+olog.dk_fname[i])
        n_dframes =  (size(dcube))[3]
        if(median_sub eq 1) then for k = 0,n_dframes-1 do dcube[*, *, k] -=  median(dcube[*, *, k])
        if(median_sub eq 2) then for k = 0, n_dframes-1 do dcube[*, *, k] -=  mean(dcube[*, *, k]) 
        ft_dcube = complexarr(dimx,dimy,n_dframes)
    endif
    print, 'Doing Fourier Transforms for file: ',filename
    for k = 0, n_frames-1 do begin
        ft_cube[*,*,k] = fft(shift(cube[*,*,k], -dimx/2,  -dimy/2)*window, 1)
        if keyword_set(disp_windowing) then begin ;;ACC: useful for checking windowing (e.g. make sure a wide companion is suppressed)
           !p.multi=[0,2,1]
           image_cont,shift(shift(cube[*,*,k], -dimx/2,  -dimy/2)*window,128,128),/n,/a,tit='Windowed'
           image_cont,shift(shift(cube[*,*,k], -dimx/2,  -dimy/2),128,128),/n,/a,tit='No Window'
           wait,0.5
        endif
        if (showps eq 1) then begin
           temp =  modsq(ft_cube[0, 0, k])/10.
           image_cont, shift(modsq(ft_cube[*,*,k]) < temp,dimx/2,dimy/2)^0.5, /noc
           wait,  0.2
        endif
     endfor

    if (olog.dk_fname[i] ne '') then for k =  0, n_dframes -1 do $
      ft_dcube[*,*,k] = fft(shift(dcube[*,*,k], -dimx/2,  -dimy/2)*window, 1)
    if (n_frames eq 1) then ps =  modsq(ft_cube) else ps = total(modsq(ft_cube),3)/float(n_frames)
    if (olog.dk_fname[i] ne '') then begin
        dps = total(modsq(ft_dcube),3)/float(n_dframes) 
        dps[0, 0] =  0.0        ;
        ft_cube[0, 0, *] -= total(ft_dcube[0, 0, *])/float(n_dframes) ;As this is used for bias subtraction in bispect
    endif else dps = fltarr(dimx,dimy)
    print, 'Now calculating bispectrum...'
    bispect, ft_cube, mf_pvct, mf_gvct, mf_ix, mf_rmat, mf_imat, v2, v2_cov, bs,bs_var, bs_cov, bs_v2_cov, $
      bl2h_ix, bs2bl_ix,bl2bs_ix, bscov2bs_ix, cp_var, bs_all, v2_all, cvis_all, closing_tri_pix=closing_tri_pix,$
      mfc_pvct=mfc_pvct, mfc_gvct=mfc_gvct, fluxes=fluxes, dark_ps=dps,n_blocks = n_blocks, cp_cov = cp_cov, $
      avar=avar, err_avar=err_avar, imsize = imsize, v2_arr=v2_arr, phs_v2corr = phs_v2corr,  $
      hole_phs = hole_phs,  hole_err_phs = hole_err_phs,hole_piston = hole_piston, subtract_bs_bias = subtract_bs_bias
    !p.multi = [0,3,1]
    ; !p.multi = [0,2,1]
    in = indgen(n_baselines)
    m = max(v2)
    ploterr, sqrt(u^2+v^2), v2, sqrt(v2_cov[in,in]), psym=5, yr = [0,min([m,1.0])],$
      ytitle = 'Raw V^2', xtitle = 'Baseline (wavelengths)'
    plothist, atan(bs, /phase)*180/!pi, xtitle = 'Closure Phase (degs)'
    if (keyword_set(nosave) eq 0) then begin
        if (i eq 0) then bs_names =  'bs'+olog.cube_fname[i,1]+ special+'.idlvar' $
        else   bs_names = [bs_names, 'bs'+olog.cube_fname[i,1]+ special+'.idlvar' ]
        ;;Make a couple of extra useful quantities from the bispect
        ;;output.
        v2_cor = cov2cor(v2_cov, sig=v2_sig)
        cp = atan(bs, /phase)
        cp_sig = sqrt(reform(bs_var[1,*])/modsq(bs))
        ;;Now save everything to file.
        save, u, v, v2, v2_sig, v2_cor, v2_cov, bs, cp, bs_var, cp_sig, bs_cov, bs_v2_cov,avar,err_avar,phs_v2corr, $
          mf_file, cp_cov,cp_var, bs_all, v2_all, cvis_all, ps, $
          filename = analysis_dir + bs_names[i]
        plog =  {mf_file:mf_file,  window_size:window_size,  coherent_vis:coherent_vis,  $
                 coherent_bs:coherent_bs,  median_sub:median_sub,  bs_names:bs_names, special:special}
        olog.logflag =  1
        save,  olog,  plog,  framestats,  filename = analysis_dir + cubeinfo_file
    endif
;Some diagnostic plots if there is anything suspicious about the data.
; if (plotall eq 0) then window,  1,  xsize = 700,  ysize = 550, ypos=50
    mf_tot=fltarr(dimx,dimy)
    bh=-1
    for j=0,n_baselines-1 do begin
        if ((bl2h_ix[0,j] ne bh) and (bl2h_ix[1,j] ne bh)) then $
          mf_tot[mf_pvct[mf_ix[0,j]:mf_ix[1,j]]]=mf_tot[mf_pvct[mf_ix[0,j]:mf_ix[1,j]]]+mf_gvct[mf_ix[0,j]:mf_ix[1,j]]
    endfor
    mf_tot=mf_tot+rotate(shift(mf_tot,-1,-1),2)
    mask =  ((mf_tot < 0.01)*100) < 0.995

; image_cont, shift(alog((ps > max(ps)*1e-7)*(1-mask)),dimx/2,dimy/2), /nocont, /asp
                                ;loadct,  0
    image_cont, shift(alog((ps > max(ps)*1e-7)*(1-mask)),dimx/2,dimy/2), /nocont, /asp 
 ;   !p.multi = 0
  ;  im = (shift(alog((ps > max(ps)*1e-7)*(1-mask)),dimx/4,dimy/4))[0:dimx/2-1,0:dimy/2-1]
  ;  image_cont, im, /nocont, /asp 

    all_ps[*,*,i] = shift(alog((ps > max(ps)*1e-7)*(1)),dimx/2,dimy/2)
;stop

if (plotall ne 0) then begin
    device,  /close
    plotfilename = analysis_dir + bs_names[i] +'.ps'
    spawn, '\mv idl.ps ' + plotfilename
    set_plot,  'x'
endif
wait,2
endfor

end

