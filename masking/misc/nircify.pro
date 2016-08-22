;This should 'nircify' all of the paths.
pro nircify, default_path

data_path="~/nirc/data/"
root_dir="~/nirc/code/"

!path=root_dir+"Speckle:"+root_dir+"Speckle/Closure_Phase_Work:"
!path=root_dir+"Misc:"+!path
!path=!path+root_dir+"Mappit/ps_anal:"
!path=!path+root_dir+"Mappit/vis_anal:"
!path=!path+root_dir+"Mappit/util:"
!path=!path+root_dir+"Imaging:"
!path=!path+root_dir+"LWS:"
!path=!path+root_dir+"LWS/Analysis:"

!path=!path+root_dir+"Jhulib:"
!path=!path+root_dir+"PTools:"
!path=!path+expand_path('+'+root_dir+'Astrolib/pro:')+":"
!path=!path+expand_path('+'+root_dir+'Deutsch:')+":"
!path=!path+root_dir+"Freudenreich:"
!path=!path+root_dir+"Minpack:"
!path=!path+default_path

;
; Enter Any Constants Here.
;
deg_pixel=0.02057*!pi/(180.0*3600.0)
sfu_pixel=1./(256.*deg_pixel)

end
