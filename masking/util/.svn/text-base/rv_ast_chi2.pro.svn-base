;; -----------------------------------------------------------------------
;; chi2 evaluation function for mixed data: RV and Masking "Astrometry"
;; questions to frantz@astro.cornell.edu
;; -----------------------------------------------------------------------

function rv_ast_chi2, params, data=data

model = data
nbAST = n_elements(model.ast)
nu = fltarr(nbAST)
rr = fltarr(nbAST)

axis   = params[0]
epsi   = params[1]
W1     = params[2]
w0     = params[3]
incl   = params[4]
t_peri = params[5]
P      = params[6]
K      = params[7]
offset = params[8]

;; ----------------------------------------------------
;; 1) calculate the residuals for the astrometry
;; ----------------------------------------------------
res1 = 0.0D

for i = 0, nbAST-1 do begin
  ;; recursive algorithm to get the excentric anomaly
  m = 2.0 * !dpi * (data.ast[i].jd - t_peri) / P
  ee = m
  for j = 0, 2000 do ee = m + epsi * sin(ee)

  ;; orbital radius and real anomaly
  rr[i] = axis * (1 - epsi * cos(ee))
  nu[i] = acos((cos(ee)-epsi)/(1-epsi*cos(ee)))

  ;; --- insure that the angle: 0 < theta <  2* pi ---
  while (ee gt 2.0 * !dpi) do ee = ee - 2.0*!dpi
  while (ee lt 0.0) do ee = ee + 2.0*!dpi

  if (ee gt !dpi) then nu[i] = -nu[i]
endfor

model.ast.asc  = rr * (cos(nu+w0) * sin(W1) + $
                       sin(nu+w0) * cos(incl) * cos(W1))

model.ast.dec = rr * (cos(nu+w0) * cos(W1) - $
                      sin(nu+w0) * cos(incl) * sin(W1))

res1 = [(data.ast.dec-model.ast.dec)/model.ast.dec_err, $
        (data.ast.asc-model.ast.asc)/model.ast.asc_err]

;; ----------------------------------------------------
;; 2 ) calculate the residuals for the radial velocity
;; ----------------------------------------------------
res2 = 0.0D
nbRV = n_elements(model.rv)

for i = 0, nbRV-1 do begin
  ;; -----------------------------------------------------------------
  ;; recursive algorithm to get the excentric anomaly th (i.e. theta)
  ;; -----------------------------------------------------------------
  m = 2.0 * !dpi * (data.rv[i].JD - T_peri) / P
  th = m
  for j = 0, 2000 do th = m + epsi * sin(th)

  ;; -----------------------------------------------------------------
  ;; now that we have theta, we can calculate the radial velocity
  ;; -----------------------------------------------------------------  

  ;; first, the 2 (x,y) components of the velocity in the orbital plane
  ;;vx = -2.0*!dpi*sin(th)/(P*(1-epsi*cos(th)))
  ;;vy =  2.0*!dpi*cos(th)*sqrt(1-epsi^2)/(P*(1-epsi*cos(th)))

  vx = -sin(th)/((1-epsi*cos(th)))
  vy =  cos(th)*sqrt(1-epsi^2)/((1-epsi*cos(th)))

  ;; then, projection on the line of sight: i.e. the radial velocity
  vxo = (vx*cos(w0)-vy*sin(w0)) * sin(incl)
  vyo = (vx*sin(w0)+vy*cos(w0)) * sin(incl)
  model.rv[i].rv = K * vyo + offset
endfor

res2 = [(data.rv.rv-model.rv.rv)/model.rv.er]

;res1[*] = 0.0
return, [res1,res2]

end
