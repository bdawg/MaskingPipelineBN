; JDM 2001Dec01		Returns Vis and Phase of Binary Disk for u,v
;           
; INPUTS:
; Model of a binary star, each with UD sizes and a ratio.
; Params:
;   params(0) = Separation (mas) 
;   params(1) = Position Angle (degs, 
;               E of N, pointing from primary -> secondary)
;   params(2) = Brightness Ratio of Primary over Secondary
;   params(3) = UD size of primary (mas)
;   params(4) = UD size of secondary (mas)
;
; u,v : units of Number of wavelengths (i.e., rad-1)
;
; outputs: 
;   visib  -- visibility
;   phases -- phases in degrees (-180,180)


pro binary_disks,u,v,params,visib,phases

delta_dec = mas2rad(params(0)*sin( (params(1)+90)*!dpi/180.))
delta_ra = -1*mas2rad(params(0)*cos( (params(1)+90)*!dpi/180.))
 ;recall +ra points at pa 90

x=sqrt(u*u + v*v)
secondary_flux = 1.0/ (1+params(2))
primary_flux = 1.0 - secondary_flux

; Primary
diameter=mas2rad(abs(params(3)))
intercept=primary_flux
index=where(x eq 0,count)
if (count gt 0) then x(index)=1e-15

      f=2*beselj(!dpi*diameter*x,1) / (!dpi*diameter*x)
if (params(3) ge 0) then f=abs(intercept*f) else f=abs(intercept/f)
vis_primary = f

; Secondary
diameter=mas2rad(abs(params(4)))
intercept=secondary_flux

      f=2*beselj(!dpi*diameter*x,1) / (!dpi*diameter*x)
if (params(4) ge 0) then f=abs(intercept*f) else f=abs(intercept/f)
vis_secondary=f

; Using boden notation (I hope!) from michelson book
phase_factor = dcomplex( cos(-2.0*!dpi*(u*delta_ra+v*delta_dec))  ,$
                         sin(-2.0*!dpi*(u*delta_ra+v*delta_dec)))

complex_vis = vis_primary + vis_secondary*phase_factor

ri2at, double(complex_vis), imaginary(complex_vis), visib, phases
phases= mod360(phases)

;stop
end
