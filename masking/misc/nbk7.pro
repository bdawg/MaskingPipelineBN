;This function will return the refractive index of BK7

function nbk7,  lambda

b =  [1.03961212, 2.317923e-1, 1.01046945]
c =  [6.00069867e-3,  2.00179144e-2,  1.03560653e2]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
