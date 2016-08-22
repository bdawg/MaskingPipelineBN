;+
; cp_editor - intended to examine scatter among a directory full of oifits files
;
; Finds all oifits in directory, then discard 'vis.oifits' and 'mrg.oifits'
;  For each file restores the closure phase, and computes robust mean and sigma 
;  Cut implemented for cubes with high clp sigma (governed by discard_sigma variable)
;  Plots clp sigma and allows interactive 'click-to-kill' 
;  Finally, all rejected data cubes are discarded and a new merged oifits file 
;   created using JDM merge_oidata routine
;  
;  median_cut keyword added. If set, then discard_sigma is interpreted as a multiplier so that
;    clp scatter larger than median*discard_sigma is rejected
;
; PGT   May 2010
; PGT   Apr 2011. Major mods. Use similar 

pro cp_editor,outfilename=outfilename,meancut_factor=meancut_factor,data_dir=data_dir,manual=manual

if(keyword_set(outfilename) eq 0) then outfilename='editmrg.oifits'
if(keyword_set(meancut_factor) eq 0) then meancut_factor=2.  ; cut files with meansquare clp this factor times the mean
if(keyword_set(data_dir) eq 0) then data_dir='./'

prefix=''
extn='.oifits'

; Find all oifits data files
file=findfile(data_dir+prefix+'*'+extn)

if file[0] eq "" then begin
    print, "No files found"
    goto, last
endif

; get rid of existing vis.oifits and mrg.oifits
gd=where( strpos(file,'vis.oifits') eq -1)
file=file(gd)
gd=where( strpos(file,'mrg.oifits') eq -1)
file=file(gd)

n_files=n_elements(file)

; Firstly accumulate vectors of mean squared clp and
; chi2 clp for data in each file

;clp_stats=fltarr(2,n_files)
meansquare_clp=fltarr(n_files)
chisquare_clp=fltarr(n_files)
filenumstr=strarr(n_files)

for i=0,n_files-1 do begin
   filename=file(i)
   extract_t3data, t3, file=filename
   ;sd=stdev(t3.t3phi,m)
   ;clp_stats[*,i] = [m,sd]
   meansquare_clp[i] = mean(t3.t3phi^2)
   chisquare_clp[i]= mean((t3.t3phi/t3.t3phierr)^2)
   filenumstr[i]=strmid(filename,strpos(filename,'.oifits')-4,4)
endfor

; First cut is made on chi2. We need a way to select an appropriate chi2
; to form a cut to reject "outliers". The method (copied from Tom Evans code)
; is actually a bit complicated and seems a bit arbitrary. Find the largest chi2 
; among that half of the data with clp's closest to zero.

smallhalf_clp = where(meansquare_clp lt median(meansquare_clp))
chi2_cut =  max(chisquare_clp[smallhalf_clp])

chi2_cull = where(chisquare_clp ge chi2_cut)


; OK now second simple cut on cases of high meansquare clp (not chi2)
meansquare_cut = median(meansquare_clp) * meancut_factor
meansq_cull = where(meansquare_clp gt meansquare_cut)

; Merge them
dummy=fltarr(n_files)
dummy[chi2_cull] =  1
if(meansq_cull[0] ne -1) then dummy[meansq_cull] = 1
cull=where(dummy gt 0)


ix=indgen(n_files)+1


; Plot-and-reject loop   
c1again:
 plot,ix,sqrt(meansquare_clp),psym=3,symsize=3
  xyouts,ix,0,filenumstr,orientation=90

  sd_clp=sqrt(meansquare_clp)

  oplot,ix[chi2_cull],sd_clp[chi2_cull],psym=1
 if(n_elements(meansq_cull) gt 1) then $
  oplot,ix[meansq_cull],sd_clp[meansq_cull],psym=7
 if(keyword_set(man_cull)) then $
  oplot,ix[man_cull],sd_clp[man_cull],psym=6

  if(keyword_set(manual) eq 0) then goto,donec1
  print,'Click on point to be removed (left of axis to finish)'
  cursor,x,y & wait,.2
  if (x lt 0) then goto,donec1

  pdist= ((ix - x)/max(ix) )^2  + $
         ((sd_clp - y)/max(sd_clp) )^2

  if(cull[0] ne -1) then pdist[cull]=1000.  ; ensure we do not find the same points twice
  w=where(pdist eq min(pdist))
  if(keyword_set(man_cull) eq 0) then man_cull=w else man_cull=[man_cull,w]
  cull=[cull,w]
goto,c1again

donec1:

if(n_elements(cull) gt 1) then ix(cull) = ix(cull) * (-1)
good=where(ix gt 0) 
; Now form a merged oidata file
merge_oidata, outfile=outfilename, infiles=file(good)

last:
stop

end

