; $Id:
; $Log: 

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
function inquire_trecs, paramname, infostr, ix

maskname =  infostr.mask
filter   =  ''

case strmid(infostr.filter[ix], 0, 3) of
  'Si1': filter =  '7.7um'
  'Si3': filter =  '9.7um'
  'Si5': filter =  '11.7um'
  'Si6':filter =  '12.3um'
  'Ope':if (strmid(infostr.optic_cfg[ix], 0, 3) eq 'Qa-') then filter =  '18.3um' else stop
  else: stop
endcase

case paramname of
  'template': return,  'trecs/mf_' + maskname + '_' + filter + '128.idlvar'
  'clog': begin
      clog =  make_clog(infostr)       
    ;  if (strcmp(maskname,  'pattern6', 8) eq 0) then clog.apply_phscorr = 0
      return,  clog
    end
  'n_blocks': return,  10
endcase

return,0

end
 

