; $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $
; $Log: inquire_lws.pro,v $
; Revision 1.4  2006/06/16 20:19:18  mireland
; Changed lots of small things. I think that most of the nirc2 changes were Peter's
; from a couple of months ago that he didn't commit. Now there is a n_blocks option
; for calc_bispect.pro (also needed in inquire) that first splits the data into
; n_blocks blocks of data before calculating variances (more realistic errors).
;
; Revision 1.3  2006/03/14 21:54:54  mireland
; Added a cal4src to clog (so that it can be modified without destroying the
; original). Added a save of quan to the final cubeinfo file.
;
; Revision 1.2  2006/03/02 21:18:05  mireland
; Fixed some bugs, and added extra sophistication to multiple calibrator
; files. Also added new LWS filters to freud etc and fixed PA.
;
; Revision 1.1  2006/01/12 23:50:23  mireland
; Commit of the new .pro calc_bispect and calibrate scripts, and the
; first LWS code version.
;
; Revision 1.2  2005/12/20 21:51:55  mireland
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and $Log: inquire_lws.pro,v $
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Revision 1.4  2006/06/16 20:19:18  mireland
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Changed lots of small things. I think that most of the nirc2 changes were Peter's
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and from a couple of months ago that he didn't commit. Now there is a n_blocks option
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and for calc_bispect.pro (also needed in inquire) that first splits the data into
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and n_blocks blocks of data before calculating variances (more realistic errors).
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Revision 1.3  2006/03/14 21:54:54  mireland
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Added a cal4src to clog (so that it can be modified without destroying the
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and original). Added a save of quan to the final cubeinfo file.
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Revision 1.2  2006/03/02 21:18:05  mireland
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Fixed some bugs, and added extra sophistication to multiple calibrator
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and files. Also added new LWS filters to freud etc and fixed PA.
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Revision 1.1  2006/01/12 23:50:23  mireland
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and Commit of the new .pro calc_bispect and calibrate scripts, and the
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and first LWS code version.
; Added $Id: inquire_lws.pro,v 1.4 2006/06/16 20:19:18 mireland Exp $ and to important files, and added the g18 mask
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
function inquire_lws, paramname, infostr, ix

maskname =  infostr.mask
filter   =  ''

case infostr.filter[0] of
  '9.9': filter =  '9.9um'
  '8.0': filter =  '8.0um'
  '12.5':filter =  '12.5um'
  '10.7':filter =  '10.7um'
  '18.75':filter =  '18.7um'
  else: 
endcase

case paramname of
  'template': return,  'lws/mf_' + maskname + '_' + filter + '.idlvar'
  'clog': begin
      clog =  make_clog(infostr)        ;For LWS, the defaults are NOT fine !!!
      if (strcmp(maskname,  'pattern6', 8) eq 0) then clog.apply_phscorr = 0
      return,  clog
    end
  'n_blocks': return,  10
endcase

return,0

end
 

