;+
; NAME:
;       READSTARLIST
;
; PURPOSE:
;       Read Keck format starlists into structure array, defined below.
;	Structure includes raw RA(hr,min,sec) Dec(sign,deg,min,sec) as
;	well as decimal degrees precessed to J2000.0
;
; CALLING SEQUENCE:
;       Result = READSTARLIST( Filename, /LGS)
;
; INPUTS:
;       Filename = Full or relative path to a Keck starlist
;
; OPTIONAL INPUT PARAMETERS:
;	/LGS - If flag is set, will return a structure array containing only
;		LGSAO targets (those with 'lgs=1' in the comment field),
;		with associated coordinates and magnitude of the first
;		tip-tilt ref. of each (any non-'lgs=1' object immediately
;		after the LGSAO target).
;
;		The Keck LGSAO format must be followed for this to work.
;		  http://www2.keck.hawaii.edu/optics/lgsao/lgsstarlists.html
;
; OUTPUTS:
;       Result - structure array as defined by variables 'ngs_str' or 
;		'lgs_str' below, depending on where /LGS flag is set.
;
; PROCEDURES USED:
;	Procedures:	PRECESS, COPY_STRUCT (IDL Astronomy User's Library)
;
; MODIFICATION HISTORY:
;       Written February 2004, A.H. Bouchez, W.M. Keck Observatory.
;	Documented and released 08/04, AHB
;	Added vmag field, /LGS flag & TTref fields. 01/05, AHB
;	Fixed bug on line 141. 01/05, AHB.
;-
FUNCTION	READSTARLIST, filename, lgs=lgs,  sval = sval

  max_char = 120		; Max length of a line
  ngs_str = { targ:	'', $	; Target name
              raj2000:	0D, $	; RA J2000, degrees
              dej2000:	0D, $	; Dec J2000, degrees
              rah:	0, $	; RA hours
              ram:	0, $	; RA minutes
              ras:	0., $	; RA seconds
              design:	'', $	; Dec sign (+1 or -1)
              ded:	0, $	; Dec degrees
              dem:	0, $	; Dec minutes
              des:	0., $	; Dec seconds
              equ:	0., $	; Equinox
              vmag:	0., $	; V magnitude
              comment:	''}	; All data after equinox

  lgs_str = { targ:	'', $	; Target name
              raj2000:	0D, $	; Target RA J2000, degrees
              dej2000:	0D, $	; Target Dec J2000, degrees
              rah:	0, $	; Target RA hours
              ram:	0, $	; Target RA minutes
              ras:	0., $	; Target RA seconds
              design:	'', $	; Target Dec sign (+1 or -1)
              ded:	0, $	; Target Dec degrees
              dem:	0, $	; Target Dec minutes
              des:	0., $	; Target Dec seconds
              equ:	0., $	; Target Equinox
              vmag:	0., $	; Target V magnitude
              comment:	'', $	; All Target data after equinox
              ttref:	'', $	; TTref name
              ttraj2000:0D, $	; TTref RA J2000, degrees
              ttdej2000:0D, $	; TTref Dec J2000, degrees
              ttrmag:	0., $	; TTref R magnitude
              ttcomment:''}	; All TTref data after equinox

;;; 1. Read file into a vector of strings

  GET_LUN,unit
  OPENR,unit,filename
  l = ''
  str = ''
  while EOF(unit) eq 0 do begin
    READF,unit,l
    str = [str,l]
  endwhile
  FREE_LUN,unit
  nl = N_ELEMENTS(str)-1
  str = str[1:nl]

;;; 2. Parse strings, recording data in structure array 'tmp'

  tmp = REPLICATE(ngs_str, nl)
  for n=0,nl-1 do begin
    tmp[n].targ = STRMID(str[n],0,16)		; first 16 char are target name
    l = STRTRIM(STRMID(str[n],16,max_char),2)	; chop off xtra spaces

    spc = STRPOS(l,' ')				; extract RA hours
    tmp[n].rah = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    spc = STRPOS(l,' ')				; extract RA minutes
    tmp[n].ram = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    spc = STRPOS(l,' ')				; extract RA seconds
    tmp[n].ras = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    tmp[n].design = '+'				; extract sign of Dec
    if STRMID(l,0,1) eq '-' then begin
      tmp[n].design = '-'			; if negative, chop off '-'
      l = STRTRIM(STRMID(l,1,max_char),2)
    endif

    spc = STRPOS(l,' ')				; extract Dec degrees
    tmp[n].ded = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    spc = STRPOS(l,' ')				; extract Dec minutes
    tmp[n].dem = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    spc = STRPOS(l,' ')				; extract Dec seconds
    tmp[n].des = STRMID(l,0,spc)
    l = STRTRIM(STRMID(l,spc,max_char),2)

    spc = STRPOS(l,' ')				; extract equinox and comments
    if spc ne -1 then begin
      tmp[n].equ = STRMID(l,0,spc)
      tmp[n].comment = STRTRIM(STRMID(l,spc,max_char),2)
    endif else tmp[n].equ = STRMID(l,0,max_char)

    c0 = STRPOS(tmp[n].comment,'vmag=') + STRLEN('vmag=')	; extract vmag
    if c0 gt (STRLEN('vmag=')-1) then begin
      c1 = STRPOS(tmp[n].comment+' ',' ',c0)
      tmp[n].vmag = FLOAT(STRMID(tmp[n].comment,c0,c1-c0))
    endif

    ra = (tmp[n].rah + tmp[n].ram/60D + $	; sum RA to decimal degrees
          tmp[n].ras/3600D ) * 15
    if tmp[n].design eq '+' then design = 1 else design = -1
    de = (tmp[n].ded + tmp[n].dem/60D + $	; sum Dec to decimal detrees
          tmp[n].des/3600D ) * design
    if tmp[n].equ ne 2000. then $		; precess to J2000 if necessary
      PRECESS, ra, de, tmp[n].equ, 2000.0
    tmp[n].raj2000 = ra
    tmp[n].dej2000 = de
  endfor

;;; 3. If /LGS flag is set, identify LGSAO targets and associate TTref data.

  if KEYWORD_SET(lgs) then begin
    slist = REPLICATE(lgs_str, nl)
    nlgs = 0
    for n=0,nl-1 do begin
      lgsflag = (STRPOS(tmp[n].comment,'lgs=1') ne -1)
      if lgsflag then begin
        COPY_STRUCT, tmp[n], lgs_str
        slist[nlgs] = lgs_str
        if n lt (nl-1) then begin		; if not last line in list...
          ttrflag = (STRPOS(tmp[n+1].comment,'lgs=1') eq -1)
          if ttrflag then begin
            slist[nlgs].ttref = tmp[n+1].targ
            slist[nlgs].ttraj2000 = tmp[n+1].raj2000
            slist[nlgs].ttdej2000 = tmp[n+1].dej2000
            slist[nlgs].ttcomment = tmp[n+1].comment
            c0 = STRPOS(slist[nlgs].ttcomment,'rmag=') + STRLEN('rmag=')
            if c0 gt (STRLEN('rmag=')-1) then begin
              c1 = STRPOS(slist[nlgs].ttcomment+' ',' ',c0)
              slist[nlgs].ttrmag = $
                FLOAT(STRMID(slist[nlgs].ttcomment,c0,c1-c0))
            endif
          endif
        endif
        nlgs = nlgs + 1
      endif
    endfor
    slist = slist[0:nlgs-1]
  endif else slist = tmp

RETURN,slist
END
