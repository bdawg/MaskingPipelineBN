;This function is for fitting to oifits data using mpffit.
FUNCTION binary_oifits_apriori, p, t3data=t3data, apriori=apriori
;  delta = 1D-2
  pfull =  double([p, 0.1, 0.1])
  modelt3 = binary_t3data(pfull,t3data=t3data)
;  q =  pfull
;  q[0] += delta
;  modelt3_0 = binary_t3data(q, t3data = t3data)
;  q =  pfull
;  q[1] += delta
;  modelt3_1 = binary_t3data(q, t3data = t3data)
;  q =  pfull
 ; q[2] += delta
 ; modelt3_2 = binary_t3data(q, t3data = t3data)
;  dp0 = [(modelt3.t3phi-modelt3_0.t3phi)/t3data.t3phierr/delta, 1./apriori.err[0], 0, 0]
;  dp1 = [(modelt3.t3phi-modelt3_1.t3phi)/t3data.t3phierr/delta, 0, 1./apriori.err[1], 0]
;  dp2 = [(modelt3.t3phi-modelt3_2.t3phi)/t3data.t3phierr/delta, 0, 0, 1./apriori.err[2]]
;  dp = [[dp0], [dp1], [dp2]]
  return,  [double(t3data.t3phi-modelt3.t3phi)/t3data.t3phierr,double(p-apriori.val)/apriori.err]
END
