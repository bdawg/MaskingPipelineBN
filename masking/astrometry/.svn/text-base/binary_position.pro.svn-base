;This function is based on Christian Hummel's 2001 Michelson School
;web document. It's designed to be used directly for fitting...
;Updates from original: vectorized multiple times
;params are... 
; T0: The time of periastron passage
; P:  The period
; a:  The semi-major axis
; e:  The eccentricity
; n:  Capital omega ( an orbital angle )
; w:  Little omega
; i:  The inclination
;NB for differentiation, see mathematica/binary_diff.nb
;Done some timing tests on my laptop: this procedure is dominated
;by calculation, not interpretation when there are more than 30 jds.

;Meaning of negative numbers:
;Negative eccentricity: change \omege by 180 degrees and T by half a
; period.
;Negative inclination: the same as positive inclination!

function binary_position, params, jds, projfact=projfact,  deriv = deriv,  vr = vr

t = jds-params[0]
P = params[1]
a = params[2]
e = abs(params[3])
n = params[4]*!pi/180.
w = params[5]*!pi/180.
i = params[6]*!pi/180.
;Allow a negative eccentricity to have the obvious meaning.
if (params[3] lt 0) then begin
 t += P/2.
 w += !pi
endif

;The mean anomaly 
;(p,t) -> M 
;Tr1 (transformation 1)
;Sequence is (p,t,e) -> M -> bE -> nu -> alpha. We want: 
;dAdp = dnudbE.dbEdM.dMdp
;dAdt = dnudbE.dbEdM.dMdp
;dAde = dnudbE.dbEde + dnude
dMdt = -2*!pi/p ;- sign because time _since_ periastron passage.
dMdp = -2*!pi*t/p^2 ;NB Can be very large.
M = 2*!pi*(t mod p)/p
;The eccentric anomaly, M = E - esinE
;(M,e) -> (bE,e) ;Tr2
;1 = dbEdM - e dbEdM cos(bE) 
;0 = dbEde - e*dbEde*cos(bE) - sin(bE)
bE = M+e*sin(M)+e^2/2*sin(2*M)
for k = 0,4 do bE=bE+(M-bE+e*sin(bE))/(1-e*cos(bE))
;The `true anomaly'. With a pi ambiguity,
;nu = 2*atan(sqrt((1+e)/(1-e))*tan(bE/2))
;(bE,e) -> (nu,e) ;Tr3
nu=2*atan(sqrt((1+e)/(1-e))*sin(bE/2), cos(bE/2))
;Offset for nu
;(nu,w) -> alpha ;Tr4
alpha=nu+w

if (arg_present(deriv)) then begin
 dbEdM = 1./(1-e*cos(bE))
 ;dedM  = (bE-1)/sin(bE) ;!!! Actually want dbEde
 dbEde =  sin(bE)/(1-e*cos(bE))
 ;Derivatives are now for alpha (just an offset from nu)
 ;From mathematica...
 dnude  = sin(bE)/(e-1)/sqrt((1+e)/(1-e))/(e*cos(bE)-1)
 dnudbE = (e-1)*sqrt((1+e)/(1-e))/(e*cos(bE) - 1)

 dAdp = dnudbE*dbEdM*dMdp
 dAdt = dnudbE*dbEdM*dMdt
 dAde = dnudbE*dbEde + dnude
endif

;Final calculations (square brackets mean trivial):
;(alpha,e,i) [+a,n] -> (rho,theta) Tr5
;We have dAd(p,t,e,w), with alpha=A. Only complex for
;e, where we need:
;drhode = drhodA.dAde + drhodnu.dAde + drhode
;dthde = dthdA.dAde + dthde
;Also, drhodp = drhodA.dAdp etc...
;!!! Much of the following can be sped-up by pre-calculation.
sqtmp = sqrt(cos(alpha)^2+sin(alpha)^2*cos(i)^2)
rho=a*(1-e^2)/(1+e*cos(nu))*sqtmp
theta=atan(sin(alpha)*cos(i),cos(alpha))+n

if (arg_present(deriv)) then begin
 drhodA = a*(e^2-1)/(1+e*cos(nu))*cos(alpha)*sin(alpha)*(sin(i))^2/sqtmp
 drhodnu = -a*e*(e^2-1)*sin(nu)/(1+e*cos(nu))^2*sqtmp
 drhode =  -a*(2*e+(1+e^2)*cos(nu))*sqtmp/(1 + e*cos(nu))^2
 drhodi = -a*(1-e^2)/(1+e*cos(nu))*cos(i)*sin(i)*(sin(alpha))^2/sqtmp
 dthdA = cos(i)/(cos(alpha))^2/(1+(cos(i)*tan(alpha))^2)
 dthdi = -sin(i)*tan(alpha)/(1+(cos(i)*tan(alpha))^2)
 ;[T0,P,a,e,n,w,i]
 zeros =  replicate(0.,  n_elements(jds))
 ones =  replicate(1.,  n_elements(jds))
 drho =  [[(drhodA+drhodnu)*dAdt], [(drhodA+drhodnu)*dAdp],  [rho/a], [drhode+(drhodA+drhodnu)*dAde], $
  [zeros], [drhodA*!pi/180.], [drhodi*!pi/180.]]
 dth  =  [[dthdA*dAdt], [dthdA*dAdp], [zeros], [dthdA*dAde], [ones*!pi/180.], [dthdA*!pi/180.], $
  [dthdi*!pi/180.]]*180/!pi
 deriv =  [[[drho]], [[dth]]]
endif

projfact = (cos(alpha))^2 ; + e*cos(w)
vr = (cos(alpha) + e*cos(w))

return, [[rho], [theta*180/!pi]]
end
