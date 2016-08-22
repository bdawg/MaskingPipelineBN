;+
; cp_scatter_cut - intended to examine scatter among a directory full of oifits files
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


pro cp_scatter_cut,outfilename,discard_sigma=discard_sigma,data_dir=data_dir,noclick=noclick,median_cut=median_cut

if(keyword_set(discard_sigma) eq 0) then discard_sigma=3.
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

clp_stats=fltarr(2,n_files)
filenumstr=strarr(n_files)

for i=0,n_files-1 do begin
   filename=file(i)
   extract_t3data, t3, file=filename
   sd=stdev(t3.t3phi,m)
   clp_stats[*,i] = [m,sd]
   
   filenumstr[i]=strmid(filename,strpos(filename,'.oifits')-4,4)
endfor


; Now edit out the ones with large CLP scatter
;  resistant_mean,clp_stats[1,*],2.0,mn0
mn0=median(clp_stats[1,*])
sd0=robust_sigma_mike(clp_stats[1,*])

cull=where(clp_stats[1,*] gt (mn0 + discard_sigma*sd0) ) 

; Here we are choosing simple cut on the median, not median+discard*sigma 
if(keyword_set(median_cut) ne 0) then $
   cull=where(clp_stats[1,*] gt (mn0 * discard_sigma) )

ix=indgen(n_files)+1

; Plot-and-reject loop   
c1again:
 plot,ix,clp_stats[1,*],psym=4
  xyouts,ix,0,filenumstr,orientation=90

 if(n_elements(cull) gt 1) then $
  oplot,ix[cull],clp_stats[1,cull],psym=2

  if(keyword_set(noclick) ne 0) then goto,donec1
  print,'Click on point to be removed (left of axis to finish)'
  cursor,x,y & wait,.2
  if (x lt 0) then goto,donec1

  pdist= ((ix - x)/max(ix) )^2  + $
         ((clp_stats[1,*] - y)/max(clp_stats[1,*]) )^2

  if(cull[0] ne -1) then pdist[cull]=1000.  ; ensure we do not find the same points twice
  w=where(pdist eq min(pdist))
  cull=[cull,w]
goto,c1again

donec1:

if(n_elements(cull) gt 1) then ix(cull) = ix(cull) * (-1)

good=where(ix gt 0) 
; Now form a merged oidata file
merge_oidata, outfile=outfilename, infiles=file(good)

last:

end

