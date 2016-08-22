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
function inquire_jwst, paramname, infostr, ix

if (infostr.mask eq '') then maskname = 'g7s6' $
 else maskname =  infostr.mask

case infostr.filter[ix] of
  'F210M': filter =  'F210N'
  'F481M': filter =  'F481M'
  'F348M': filter =  'F348M'
  'F403M': filter =  'F403M'
  else: filter = infostr.filter[ix]
endcase 

case paramname of
  'template': begin
     datasize=fix(infostr.cube_sz[ix,0])
     if(datasize eq 256) then szname='' else szname='_'+strtrim(string(datasize),2)
      return,  'jwst/mf_' + maskname + '_' + filter + szname + '.idlvar'
    end
  'clog': begin
      clog =  make_clog(infostr)       
    ;  if (strcmp(maskname,  'pattern6', 8) eq 0) then clog.apply_phscorr = 0
      return,  clog
    end
  'n_blocks': return,  10
endcase

return,0

end
 

