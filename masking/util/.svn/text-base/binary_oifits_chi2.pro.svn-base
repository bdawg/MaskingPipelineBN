;This function is for fitting to oifits data using mpffit.
;!!! Need to change binary_vis2data to
;binary_vis2atm... Unfortunately, the hessian matrix etc doesn't
;really work for binary_vis2atm and needs more thought !!! (c.f. JWST emails)

FUNCTION binary_oifits_chi2, p,  hessian = hessian,  covar = covar,  corr = corr,  sig = sig, print=print
common t3block,  t3data,  apriori, vis2data, usevis,  cp_cinv,  proj,  proj_err, force_min_sep
 
  delpar =  [0.1, 0.1, 0.01*p[2]]
  pfull =  [p, 0.1, 0.1]
  if (usevis) then modelv2 = binary_vis2data(double(pfull),vis2data=vis2data)
  modelt3 = binary_t3data(double(pfull),t3data=t3data)
  t3resid = double(mod360(t3data.t3phi-modelt3.t3phi))
  if (proj[0] ne -1) then begin
   chi2 = total((t3resid#proj)^2/proj_err^2)
  endif else if (cp_cinv[0] eq -1) then chi2 = total([t3resid/t3data.t3phierr]^2) $
  else chi2 = transpose(t3resid)#cp_cinv#t3resid
  if keyword_set(print) then begin
   print, chi2
  endif
  chi2 +=  total([double(p[0:1]-apriori.val[0:1])/apriori.err[0:1]]^2)
  chi2 += (  double(1./p[2]-1./apriori.val[2])*double(apriori.val[2])^2./apriori.err[2]  )^2
  if (p[0] lt force_min_sep) then chi2 += 1e3
  if (force_min_sep gt 0 and p[2] lt 3) then chi2 += 1e3
  if (usevis) then chi2 += total([double(vis2data.vis2data-modelv2.vis2data)/vis2data.vis2err]^2)
  if (arg_present(hessian) or arg_present(covar) or arg_present(sig)) then begin
    nt3 = n_elements(t3data)
    nv2 = n_elements(vis2data)
    ;;We are looking for second derivatives of chi^2
    ;;with respect to the model parameters (sep, pa, contrast)
    ;;2*dp(p-apriori.val)/apriori_err is
    ;;the apriori component (diagonal). The other component is
    ;;total(2/t3phierr^2*Dt3phi_iDt3phi_j) 

    ;;jacob is the derivative of V^2 or closure-phase with respect to 
    ;;the three parameters.
    jacob =  dblarr(nt3, 3)
    jacob_v = dblarr(nv2, 3)
    for i = 0, 2 do begin
     q = pfull
     q[i] += delpar[i]
     modelt3_1 = binary_t3data(double(q),t3data=t3data)
     jacob[*, i] = mod360(modelt3_1.t3phi -modelt3.t3phi)/delpar[i]
     if (usevis) then begin 
      modelv2_1 = binary_vis2data(double(q),vis2data=vis2data)
      jacob_v[*, i] = (modelv2_1.vis2data -modelv2.vis2data)/delpar[i]
     endif
    endfor
    hessian = diag_matrix(2./apriori.err^2)
    for i = 0, 2 do for j = 0, 2 do begin
      if (proj[0] ne -1) then begin
       hessian[i, j] += total((jacob[*, i]#proj)*(jacob[*, j]#proj)*2/proj_err^2)
      endif else if (cp_cinv[0] eq -1) then $
       hessian[i, j] += total(jacob[*, i]*jacob[*, j]*2/t3data.t3phierr^2) $
      else hessian[i, j] += transpose(jacob[*, i])#cp_cinv#jacob[*, j]*2
      
      if (usevis) then hessian[i, j] += total(jacob_v[*, i]*jacob_v[*, j]*2/vis2data.vis2err^2)
    endfor
    covar =  invert(hessian)*2.
    corr = cov2cor(covar,  sig = sig)
  endif
  return,  chi2
END
