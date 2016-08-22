;This procedure plots the phase of a binary system, with the phase
;extracted from closure phase over-plotted
;...test parameters...
;For GJ 84, try p=[420,47,15,0.1,0.1] !!! Not so good...
;For GJ 164 try p=[88.4, 100.0 ,5.13, 0.1, 0.1]

pro plot_binary_phase, cubeinfo_file,  root_dir, p

restore,  cubeinfo_file
if (olog.logflag ne  2) then begin
 print,  'Error: Data hasn''t been calibrated yet!'
 return
end
binary_disks,u1,v1,p,visib,phases
mf_filestring = root_dir + '/templates/' + plog.mf_file
restore, mf_filestring
w = where(quan.phase_err gt 0, complement=w2)
known_phase = w2
order_known_phase, bl2h_ix, known_phase, known_dir

;Make some phasors...
mvis = exp(complex(0,phases*!pi/180))
vis = exp(complex(0, quan.phase))
hole_phasors = complexarr(n_holes)
hole_phasors[*] = 1.0
for i = 0,(n_elements(known_phase) - 1) do begin
 if (known_dir[i] eq 0) then $;forwards
   hole_phasors[bl2h_ix[1,known_phase[i]]] = hole_phasors[bl2h_ix[0,known_phase[i]]]*mvis[known_phase[i]] $
 else $ ;backwards
   hole_phasors[bl2h_ix[0,known_phase[i]]] = hole_phasors[bl2h_ix[1,known_phase[i]]]*conj(mvis[known_phase[i]])
endfor
new_vis = vis*conj(hole_phasors[bl2h_ix[0,*]])*hole_phasors[bl2h_ix[1,*]]
new_phase =  atan(new_vis,  /phase)

vp = cos(p[1]*!pi/180.)
up = sin(p[1]*!pi/180.)
m =  max(sqrt(u1^2+v1^2))
uplot =  (findgen(200)/100. - 1.0)*up*m
vplot =  (findgen(200)/100. - 1.0)*vp*m
binary_disks,uplot,vplot,p,dummy,plot_phases
m2 =  max(abs(plot_phases))
plot, uplot*up+vplot*vp, plot_phases,  xstyle = 1,  yr = [-1.5*m2,  1.5*m2],  xtitle = 'Spatial Frequency',  ytitle = 'Phase (degs)'
oploterr, u1[w]*up+v1[w]*vp, new_phase[w]*180/!pi, quan.phase_err[w]*180/!pi, psym=4
oplot,  u1[w2]*up+v1[w2]*vp, new_phase[w2]*180/!pi, psym=5

end
