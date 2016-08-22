;$Id: checkpharospeckle.pro,v 1.1.1.1 2005/12/19 05:15:05 mireland Exp $
;$Log: checkpharospeckle.pro,v $
;Revision 1.1.1.1  2005/12/19 05:15:05  mireland
;Initial release of masking code.
;
;Revision 1.1  2003/09/14 11:02:43  jpl
;add specklecube mode
;

; pro checkpharospeckle, fn
   frames=readfitspharo(fn,header,/specklecube)
   
   object=sxpar(header,'OBJECT')
   filter=sxpar(header,'FILTER')
   lyot  =sxpar(header,'LYOT')
   carousel=sxpar(header,'CAROUSEL')

   print, fn+'    '+object+'    '+filter+' '+lyot+' '+carousel

   nframes= (size(frames))[3]

   tot = fltarr(nframes)
   mins = fltarr(nframes)
   maxs = fltarr(nframes)
   sd  = fltarr(nframes)
   cummax = fltarr(nframes)
   
   sum=0
   for i=0,nframes-1 do begin
      f = frames[*,*,i]
      tot[i] = total(f)
      mins[i] = min(f)
      maxs[i] = max(f)
      sd[i] = sqrt( (moment(f))[1] )  
      sum = sum + maxs[i]
      cummax[i] = sum
      tvscl,f
   endfor
   
   !p.multi=[0,2,2]
   plot, tot,title='total counts per frame'
   plot, maxs,title='peak per frame'
   plot, sd,title='st dev'
   plot, cummax, title='cumulative max'
   
   
 end
