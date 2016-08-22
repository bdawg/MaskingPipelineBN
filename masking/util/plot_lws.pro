;This is a special plotting program for LWS data, which comes in
;sets of 4 oiftis files that have to be concatenated.

pro plot_lws, dir,  filenum,  subarrs,  v2data = v2data,  rad = rad,  rv2 = rv2,  err_rv2 = err_rv2

outfile =  dir+filenum[0]+'all.oifits'
merge_oidata, outfile=outfile, infiles=dir+filenum+subarrs+'.oifits'
extract_vis2data, v2data, file=outfile
ave_vis,  sqrt(v2data.u^2 +v2data.v^2) ,v2data.vis2data,v2data.vis2err,rad,rv2,err_rv2,  auto = 0.2/10e-6
ploterr,  rad,  rv2,  err_rv2,  psym = 4,  yr = [0, 1.1],  ytitle = 'V^2',  xtitle = 'Baseline (wavelengths)'

end
