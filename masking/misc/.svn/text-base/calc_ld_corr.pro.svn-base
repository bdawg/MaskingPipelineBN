;This function returns the limb-darkening correction based on the
;4-parameter models of claret (2000). Based on V^2 of 0.3
;I could download table 15 from the web if I was keen...

function calc_ld_corr, a1, a2, a3, a4

mu = findgen(1000)/1000.0
f = 1 - a1*(1-mu^0.5) -a2*(1-mu) - a3*(1-mu^1.5) - a4*(1-mu^2)
;plot, sqrt(1.0-mu^2), f
v = hankel(reverse(sqrt(1.0-mu^2)), reverse(f), u=u, maxu=1.0)

corr  = (where(v lt sqrt(0.3)))[0]
plot, u, v^2
return, interpol(u,v^2,0.3)/0.661878
end
