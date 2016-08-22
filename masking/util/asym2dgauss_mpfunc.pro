;

function asym2dgauss_mpfunc, params, dp, vis2data0=vis2data0,t3data0=t3data0,$
                           vis2model0=vis2model0,t3model0=t3model0,$
				scale=scale,dim=dim, $
 			    x=x,y=y,squarex=squarex,squarey=squarey,$
 			aa=aa,tt=tt,im=im


;print,'Params: ',params
; note vis2data/t3data are from extract_vis2data calls, not the
;  original binary tables.
;vis2model0=vis2data0
vis2model0 = asym2dgauss_vis2data(params,vis2data=vis2data0,$
 			scale=scale,dim=dim, $
                        x=x,y=y,squarex=squarex,squarey=squarey,$
                        aa=aa,tt=tt)


;t3model0=t3data0

t3model0=asym2dgauss_t3data(params,t3data=t3data0,$
                        scale=scale,dim=dim, $
                        x=x,y=y,squarex=squarex,squarey=squarey,$
                        aa=aa,tt=tt,im=im)



residuals = [ (vis2data0.vis2data - vis2model0.vis2data)/(vis2data0.vis2err),$
              angle_diff(t3data0.t3phi,t3model0.t3phi)/t3data0.t3phierr]

return,residuals


end

