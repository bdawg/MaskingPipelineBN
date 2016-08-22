; ---- Usage: hlog=lacan(datadir,fitsfilename)
; Parses VAMPIRES logfiles and returns a structure of useful info
; that otherwise would come from FITS headers.
; 
; If no filename is given use 'LogfileVampires.txt'
;
; Note - this system was used in 2013 only.

function lacan,datadir,fitsfilename,logfilename=logfilename,gen2filename=gen2filename

if keyword_set(logfilename) eq 0 then logfilename='LogfileVampires.txt'
if keyword_set(gen2filename) eq 0 then gen2filename='LogfileGen2.txt'

s='' &  i=0 &  f=float(0.0) &  d=double(0.0)

headinfo={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,optic_cfg:s, $
          lyot:s,grism:s,source_name:s,utc:s,date:s,jd:d,elevation:f,airmass:f,pa:f, $
          ra:d,dec:d,equinox:f, mask:s,  raoff:f, decoff:f,  del_elev:f,  del_pa:f, $
          emgain:f, hwp:f, timingpattern:s, imgRotAng:f, imgRotPad:f, imgRotPap:f ,$
          ADCstagepos:f, ADCp1angle:f, ADCp2angle:f, azimuth:f, localtime:s} ;These last 2 lines are added for VAMPIRES

headinfo.instrument='VAMPIRES'

;print,'### Lacan skipped with dodgy hack - April 2014 testing ###'
;goto,skipeverything

; Part 1 - Read VAMPIRES logfile
fullFileName=[datadir+logfilename]
formats='A,A,A,A,F,F,A,A,A'
readcol,fullFileName,F=formats,delimiter=' ',ldate,ltime,lfitsfilepath,ltimingmode,lemgain,lhwp,loldmaskfield,lfilter,lmask

searchstring=['*'+fitsfilename]
curRow=where(strmatch(lfitsfilepath,searchstring))

headinfo.filter=lfilter[curRow]
headinfo.localtime=ltime[curRow]
curTime=ltime[curRow]
headinfo.date=ldate[curRow]
headinfo.emgain=lemgain[curRow]
headinfo.hwp=lhwp[curRow]
headinfo.mask=lmask[curRow]
headinfo.timingpattern=ltimingmode[curRow]




; Part 2 - Read Gen2 automatic logfile
fullG2FileName=[datadir+gen2filename]
openr,lun,fullG2FileName,/get_lun

; Loop through file until a matching time (to second) is found.
; Assume 3 lines, always in the order IMGROT, ADC, TELESCOPE
success=0
searchstring=['*'+curTime+'*']
while success eq 0 do begin
    line=''
    readf,lun,line
    if strmatch(line,searchstring) then success = 1
end
; Check this line is actually a SCEXLOG line (not some other event)
if strmatch(line,'*SCEXLOG*') NE 1 then begin
    line=''
    readf,lun,line
endif

; Read SCEXLOG:IMGROT line
pos=strpos(line,'angle')
pos=pos+6
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.imgRotAng=val

pos=strpos(line,'pad')
pos=pos+4
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.imgRotPad=val

pos=strpos(line,'pap')
pos=pos+4
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.imgRotPap=val


; Read SCEXLOG:ADC line
line=''
readf,lun,line

pos=strpos(line,'ra=')
pos=pos+3
valStr=strmid(line,pos,10)
headinfo.ra=ten(valStr)

pos=strpos(line,'dec=')
pos=pos+4
valStr=strmid(line,pos,10)
headinfo.dec=ten(valStr)

pos=strpos(line,'pa=')
pos=pos+3
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.pa=val

pos=strpos(line,'stagepos=')
pos=pos+9
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.ADCstagepos=val

pos=strpos(line,'p1angle=')
pos=pos+8
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.ADCp1angle=val

pos=strpos(line,'p2angle=')
pos=pos+8
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.ADCp2angle=val

; Read SCEXLOG:TELESCOPE line
line=''
readf,lun,line

pos=strpos(line,'AZ = ')
pos=pos+5
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.azimuth=val

pos=strpos(line,'EL = ')
pos=pos+5
valStr=strmid(line,pos,8)
reads,valStr,val
headinfo.elevation=val

free_lun,lun

skipeverything:
return,headinfo

end
