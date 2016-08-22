
restore,'TestingAdjTempl_testps_256_9h_ib224.idlvar'

restore,'TestingAdjTempl_mf_9Holes_IB_2.24_pgt_S13.idlvar'
msk_rot = 0.0*!pi/180
scale=1.0
rot1=[[cos(msk_rot),sin(msk_rot)],[-sin(msk_rot),cos(msk_rot)]]
xy_coords = scale*rot1##xy_coords

chipsz=(size(pspc))[1]
ps=shift(pspc,chipsz/2,chipsz/2)

; Firstly fit the gaussian for each splodge
splodge_fits=fltarr(7,n_baselines)

for bl=0,n_baselines-1 do begin
     thisholepair=bl2h_ix[*,bl]
     xspot=-(xy_coords[bl2h_ix[0,bl],0] - xy_coords[bl2h_ix[1,bl],0])/filter[0]*rad_pixel*chipsz+chipsz/2
     yspot=-(xy_coords[bl2h_ix[0,bl],1] - xy_coords[bl2h_ix[1,bl],1])/filter[0]*rad_pixel*chipsz+chipsz/2
     xspot_int=nint(xspot)
     yspot_int=nint(yspot)
     mf=fltarr(chipsz,chipsz)
     cookiecutter,mf,xspot_int,yspot_int,sampledisk_r
     thisspot_c=where(mf gt 0)                                 ; we make up 2 samplings. This 'c' one is centered at [chipsz/2,chipsz/2]

     startpk=max(ps * mf)
     guess_spotsize=hole_diam/filter[0]*rad_pixel*chipsz/1.9
  stop 
     this_splodge_fit=mpfitfun('gauss_splodge_fit_fn',mf,ps * mf,fltarr(chipsz,chipsz)+1,[startpk,guess_spotsize,0.,xspot,yspot]);,/quiet)
     splodge_fits[*,bl]=[this_splodge_fit,xspot,yspot]
endfor

splodge_fits[3:6,*] -= chipsz/2

plot,splodge_fits[3,*],splodge_fits[4,*],psym=4,/isotropic
oplot,splodge_fits[5,*],splodge_fits[6,*],psym=1

rad_measured=sqrt(splodge_fits[3,*]^2 + splodge_fits[4,*]^2)
rad_predict=sqrt(splodge_fits[5,*]^2 + splodge_fits[6,*]^2)

plate_distortion=median((rad_measured-rad_predict)/rad_predict)
print,' Plate Scale Adjustment =',plate_distortion

theta_measured=atan(splodge_fits[3,*],splodge_fits[4,*])
theta_predict=atan(splodge_fits[5,*],splodge_fits[6,*])
angle_error=median(theta_measured-theta_predict)
print,' Angle Adjustment =',angle_error,' radians or ',angle_error*180./!pi,' degrees'


end
