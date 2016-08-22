
; $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $
; $Log: inquire_pharo.pro,v $
; Revision 1.15  2007/06/14 18:59:19  mireland
; Some of the nirc2 stuff shouldn't have changed. Trying to finish the
; cp_cov stuff.
;
; Revision 1.14  2007/06/12 23:52:41  mireland
; Lots of changes: improved cosmic ray rejection, better naming and
; directory stuff, modifications to closure-phase histogram plots...
;
; Revision 1.13  2007/06/01 11:53:41  mireland
; Now inquire knows about different PHARO chip sizes.
;
; Revision 1.9  2006/06/16 20:21:55  mireland
; Only the inquire_ files should have changed. These now return a default
; 0 for n_blocks, but are untested.
;
; Revision 1.8  2006/03/14 21:54:54  mireland
; Added a cal4src to clog (so that it can be modified without destroying the
; original). Added a save of quan to the final cubeinfo file.
;
; Revision 1.7  2006/03/01 22:14:37  frantzm
; 9h, Ks matched filter added
;
; Revision 1.6  2006/01/12 23:50:23  mireland
; Commit of the new .pro calc_bispect and calibrate scripts, and the
; first LWS code version.
;
; Revision 1.2  2005/12/20 21:51:55  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and $Log: inquire_pharo.pro,v $
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.15  2007/06/14 18:59:19  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Some of the nirc2 stuff shouldn't have changed. Trying to finish the
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and cp_cov stuff.
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.14  2007/06/12 23:52:41  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Lots of changes: improved cosmic ray rejection, better naming and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and directory stuff, modifications to closure-phase histogram plots...
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.13  2007/06/01 11:53:41  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Now inquire knows about different PHARO chip sizes.
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.9  2006/06/16 20:21:55  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Only the inquire_ files should have changed. These now return a default
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and 0 for n_blocks, but are untested.
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.8  2006/03/14 21:54:54  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Added a cal4src to clog (so that it can be modified without destroying the
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and original). Added a save of quan to the final cubeinfo file.
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.7  2006/03/01 22:14:37  frantzm
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and 9h, Ks matched filter added
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Revision 1.6  2006/01/12 23:50:23  mireland
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and Commit of the new .pro calc_bispect and calibrate scripts, and the
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and first LWS code version.
; Added $Id: inquire_pharo.pro,v 1.15 2007/06/14 18:59:19 mireland Exp $ and to important files, and added the g18 mask
; for PHARO's inquire.
;
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
function inquire_pharo, paramname, infostr, ix

camera=infostr.instrument[ix]
maskname =  infostr.mask
filter   =  ''

case infostr.grism[ix] of
  'Fe II 1.643': begin
    filter =  'feii'
    lambda = 1.648
    end
  'H-cont 1.67' : begin
    filter =  'hcont'
    lambda = 1.668
    end
  else:
endcase

case infostr.filter[ix] of
  'Br-gamma': begin
    filter =  'brg' 
    lambda = 2.18e-6
    end
  'H': begin
    filter =  'h' 
    lambda = 1.635e-6
    end
  'J': begin
    filter =  'j' 
    lambda = 1.246e-6
    end
  'K_short': begin
    filter = 'ks' 
    lambda = 2.145e-6
    end
  'K': begin
    filter = 'k' 
    lambda = 2.196e-6
    end
  'CH4_S' : begin
    filter = 'ch4s' 
    lambda = 1.570e-6
    end
  else: 
endcase

;!!! Generally, this shouldn't be required, as the mask name is set
;automatically in freud. These lines should probably be deleted
case infostr.mask of
  '18': maskname =  'p18'
  'g18' :  maskname =  'p18'
  else:
endcase

case infostr.optic_cfg[ix] of
  '25 mas': pscale = 25.09 ;From Stanimir Metchev's document
  '40 mas': pscale = 39.93 ;From the PHARO manual.
  '-1': pscale = 25.09 ;For back-compatability
endcase

chipsz = fix(infostr.cube_sz[ix, 0])
if (chipsz ne 256) then special = strtrim(chipsz, 2) else special = ''
clog = make_clog(infostr)
clog.average_src = 1 ;This is the default for PHARO (equitorial)

case paramname of
  'template': return,  'pharo/mf_' + maskname + '_' + filter + special + '.idlvar'
  'clog': return, clog
  'n_blocks': return,  20
  'lambda': return,  lambda
  'pscale':return,  pscale
endcase

return,0

end
 

