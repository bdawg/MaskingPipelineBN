
;+
; inquire(paramname,infostr,ix=ix)
;
; inquire will return a default setting (for example, a template) by
;   using the observing log info structure.

; Input:  name         - the name of the thing you want
;         ix           - index for addressing infostr(if needed)
; Returns: infostr     - a structure with everything you wanted to know
; 
; created                                                    PGT 22Oct05
;
; This is really example code to show how to build it at
; present
function inquire_nirc2,paramname,infostr,ix

camera=infostr.instrument[ix]

maskname =  infostr.mask
filter   =  ''

case infostr.filter[ix] of
  'J + 9holeMsk': begin
    filter =  'j'
    mask = 'g9'
    end
 'H + 9holeMsk': begin
   filter =  'h'
    mask = 'g9'
   end
  'K + 9holeMsk': begin
    filter =  'k'
    mask = 'g9'
    end
  'Ks + 9holeMsk': begin
    filter =  'ks'
    mask = 'g9'
    end
  'Kp + 9holeMsk': begin
    filter =  'kp'
    mask = 'g9'
    end
  'Ks + 9holeMsk': begin 
    filter =  'ks'
    mask = 'g9'
    end
  'H2O_ice + 9holeMsk': begin
    filter =  'h2o'
    mask = 'g9'
    end
  'Lp + 9holeMsk': begin
    filter =  'lp'
    mask = 'g9'
    end
  'Ms + 9holeMsk': begin
    filter =  'ms'
    mask = 'g9'
    end
  '18hole_msk + Jcont': begin
    filter =  'jcont'
    mask = 'g18'
    end
  '18hole_msk + Hcont': begin
    filter =  'hcont'
    mask = 'g18'
   end
  '18hole_msk + Kcont': begin
    filter =  'kcont'
    mask = 'g18'
    end
  '18hole_msk + PaBeta': begin
    filter =  'pabeta'
    mask = 'g18'
    end
  'CH4_short + 9holeMsk': begin
    filter = 'ch4s'
    mask = 'g9'
    end
   else:stop
endcase

case paramname of
  'template': return,  'nirc2/mf_' + maskname + '_' + filter + '.idlvar'
  'clog': return,  make_clog(infostr)  ;For NIRC2, the defaults are fine.
  'n_blocks': return,  0
  'pscale' : return,  9.97 ;!!! Actually multiple scales...
  'mask':return,  mask 
 endcase
return,0

end
 

