PRO lookat_chisqmap,filename



restore,filename

;Scale chisq - assumes model should have fit with chisq=1
;chisqmap=chisqmap/min(chisqmap)

chisqmap=chisqmap*dof ;Go back to raw chisq values

;Show chisq map
loadct,39
!p.thick=1.0
!p.charthick=1.0
!p.charsize=1.4
image_cont,chisqmap,xval=a0vals,yval=a5vals,xtit='Dust shell radius (mas)',ytit='Dust-scattred intensity fraction',/noc,tit=filename;tit='W Hya 2010 1.24 microns'
;contour,chisqmap,a0vals,a5vals,/overplot,levels=[min(chisqmap)+1] ;1 sigma contour

;surface,chisqmap,a0vals,a5vals,xtit='Dust shell radius (mas)',ytit='Dustfract',charsize=3.0

;Get 1 sigma errors... this way is too complicated...
; Chisq is -log of probability density function
; a0proj=total(chisqmap,2)
; pdf_a0=exp(-a0proj)
; pdf_a0=pdf_a0/total(pdf_a0)
; cdf_a0=fltarr(n_elements(pdf_a0))
; cdf_a0[0]=pdf_a0[0]
; for i = 1,n_elements(cdf_a0)-1 do cdf_a0[i]=cdf_a0[i-1] + pdf_a0[i]
; ;+/-1sigma is where pdf is 15.87 and 84.14 
; a5proj=total(chisqmap,1)
; pdf_a5=exp(-a5proj)
; pdf_a5=pdf_a5/total(pdf_a5)
; cdf_a5=fltarr(n_elements(pdf_a5))
; cdf_a5[0]=pdf_a5[0]
; for i = 1,n_elements(cdf_a5)-1 do cdf_a5[i]=cdf_a5[i-1] + pdf_a5[i]


;Use a quintic interpolated map, to find a minimum and errors
chisqmap_interp=tri_surf(chisqmap,xvalues=a0vals,yvalues=a5vals,nx=100,ny=100)
;chisqmap_interp=congrid(chisqmap,100,100)
locate_peak,-chisqmap_interp,xpk,ypk
a0vals_interp=interpol(a0vals,100)
a5vals_interp=interpol(a5vals,100)
a0min=a0vals_interp(xpk)
a5min=a5vals_interp(ypk)

;image_cont,chisqmap_interp,xval=a0vals_interp,yval=a5vals_interp,xtit='Dust shell radius (mas)',ytit='Dustfract',/noc,tit=filename
contour,chisqmap_interp,a0vals_interp,a5vals_interp,/overplot,levels=[min(chisqmap)+1]
oplot,[a0min],[a5min],psym=1

;Find 1 sigma errors (do lower bounds and upper bounds separately)
chicut=min(chisqmap_interp)+1
chicut=min(chisqmap)+1 ;Are error bounds based on interpolated of raw data? Hopefully this effect dissaperas with zoomed-in maps
a5lb_ind=closest(chicut,chisqmap_interp[xpk,0:ypk],decide=upper)
a5ub_ind=closest(chicut,chisqmap_interp[xpk,ypk:99],decide=upper)+ypk
a0lb_ind=closest(chicut,chisqmap_interp[0:xpk,ypk],decide=upper)
a0ub_ind=closest(chicut,chisqmap_interp[xpk:99,ypk],decide=upper)+xpk
oplot,[a0vals_interp[a0lb_ind],a0vals_interp[a0ub_ind]],[a5vals_interp[ypk],a5vals_interp[ypk]],psym=2
oplot,[a0vals_interp[xpk],a0vals_interp[xpk]],[a5vals_interp[a5lb_ind],a5vals_interp[a5ub_ind]],psym=2

; Scale errors by sqrt(chisq/dof)
a0_loerr=(a0min-a0vals_interp[a0lb_ind])*sqrt(min(chisqmap_interp)/dof)
a0_hierr=(a0vals_interp[a0ub_ind]-a0min)*sqrt(min(chisqmap_interp)/dof)
a5_loerr=(a5min-a5vals_interp[a5lb_ind])*sqrt(min(chisqmap_interp)/dof)
a5_hierr=(a5vals_interp[a5ub_ind]-a5min)*sqrt(min(chisqmap_interp)/dof)

;Only report symmtric errors, easier later, close enough for now.
print,'Dust shell radius:',a0min,' +/-',(a0_hierr+a0_loerr)/2
print,'Dustfract:',a5min,' +/-',(a5_hierr+a5_loerr)/2
print,'Reduced CHISQ at minimum:',chisqmap_interp[xpk,ypk]/dof

amx=n_elements(a0vals)-1
xyouts,a0vals[0]-(a0vals[amx]-a0vals[0])/5,a5vals[0]-(a5vals[amx]-a5vals[0])/8,[string(a0min)+' +/-'+strtrim(string((a0_hierr+a0_loerr)/2),2)]
xyouts,a0vals[8]-(a0vals[amx]-a0vals[0])/7,a5vals[0]-(a5vals[amx]-a5vals[0])/8,[string(a5min)+' +/-'+strtrim(string((a5_hierr+a5_loerr)/2),2)]


end

