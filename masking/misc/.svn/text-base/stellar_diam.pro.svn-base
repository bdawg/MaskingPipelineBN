;This is a function that uses 2 colour photometry to predict the
;diameter of a star.
;
;TODO: Add an option for a spectral type, and de-redden the vmag. With
;a Vmag correction of 3.1*delta(B-V), we get an optical depth 
;proportional to lambda^(-1.25). 



function stellar_diam, mag_blue, mag_red, method=method

if(keyword_set(mag_blue) eq 0) then begin
   print,'diam = stellar_diam(mag_blue, mag_red, method=x)'
   print,'methods are:'
   print,'0 - Davis B and V mag relationship'
   print,'1 - van Belle 99 V and K mag for Evolved Stars'
   print,'2 - van Belle 99 B and K mag for Evolved Stars'
   print,'3 - van Belle 99 V and K mag for Miras'
   print,'4 - van Belle 99 B and K mag for Miras'
   print,'5 - van Belle 99 V and K mag for Main Sequence'
   print,'6 - van Belle 99 B and K mag for Main Sequence'
   return,''
endif

if(keyword_set(method) eq 0) then method=0

case method of

  0: return, (6.477*(mag_blue-mag_red) + 3.04) * 10^(-mag_red/5.0)

  1: return, 10^(.669 + .223*(mag_blue-mag_red))/10^(mag_blue/5.0)

  2: return, 10^(.648 + .220*(mag_blue-mag_red))/10^(mag_blue/5.0)

  3: return, 10^(.789 + .218*(mag_blue-mag_red))/10^(mag_blue/5.0)

  4: return, 10^(.840 + .211*(mag_blue-mag_red))/10^(mag_blue/5.0)

  5: return, 10^(.500 + .264*(mag_blue-mag_red))/10^(mag_blue/5.0)

  6: return, 10^(.500 + .290*(mag_blue-mag_red))/10^(mag_blue/5.0)

endcase 

end
