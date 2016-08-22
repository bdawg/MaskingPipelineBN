; This function defines the olog data structure for use
; in raw data analysis....
;
; Version 0.0 From cube_nirc   PGT 23Dec05
;$Id: make_olog.pro,v 1.9 2010/06/18 05:08:55 snert Exp $
;$Log: make_olog.pro,v $
;Revision 1.9  2010/06/18 05:08:55  snert
;Added quadrant
;
;Revision 1.8  2010/04/01 04:45:50  snert
;Added release clip line
;
;Revision 1.7  2009/05/28 05:34:38  snert
;Included flat_field in the olog structure.
;
;Revision 1.6  2007/06/12 23:52:41  mireland
;Lots of changes: improved cosmic ray rejection, better naming and
;directory stuff, modifications to closure-phase histogram plots...
;
;Revision 1.5  2007/06/09 22:31:43  mireland
;Lots of bug fixes and minor feature adds (mainly for using pharomkcube for full
;pupil image analysis).
;
;Revision 1.4  2007/05/30 22:43:09  frantzm
;
;added "frames_ix" in olog: begin/end frames making a datacube.
;
;Revision 1.3  2006/06/16 20:19:18  mireland
;Changed lots of small things. I think that most of the nirc2 changes were Peter's
;from a couple of months ago that he didn't commit. Now there is a n_blocks option
;for calc_bispect.pro (also needed in inquire) that first splits the data into
;n_blocks blocks of data before calculating variances (more realistic errors).
;
;Revision 1.2  2006/01/06 19:06:39  mireland
;Added uflip to olog, created a better name for cube_pharo and
;added documentation.
;
;_____________________________________________________________
;######## RELEASE CLIP HERE ######## 
function make_olog,n_cubes,tsize,frames,skies

; Define data structure to hold Header info and comments
s=replicate('',n_cubes) &  i=intarr(n_cubes) &  f=fltarr(n_cubes) &  d=dblarr(n_cubes)
fi = intarr(2, n_cubes)
olog={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,		$ ; from headers
			slit:s,optic_cfg:s,lyot:s,grism:s,source_name:s,utc:s,				$ ; ..
			date:s,jd:d,elevation:f,del_elev:f,airmass:f,pa:f,del_pa:d,		$ ; ..
      ra:d,dec:d,equinox:f,  ha:d,raoff:d,  decoff:d,               $ ; ..
      tsize:tsize,mask:'',frames:frames,skies:skies,rawdir:'',			$ ; User script input
      cubedir:'', frames_ix:fi,                                     $ ; begin/end frames of cubes
      adate:'',comments:'',proc_hist:'',logflag:0,                  $ ; ..
      cal4src:intarr(n_cubes,n_cubes), quad:i,                      $ ; ..
      dk_fname:s,cube_fname:strarr(n_cubes,2),cube_tsize:f,         $ ; output this script
      cube_sz:fltarr(n_cubes,3),  uflip:1.0, flat_file:strarr(1)}   	                                  ;

return,olog
end
