;Call this procedure to display images from MACIM.
;The idea is that this should be optimised for multi-thread (expert)
;diagnostics. The python gui is for simple single-thread stuff.
;
;NB This must be run from the same directory as MACIM directory (cd, 'DIR')
pro display_macim

window,  0
window,  1

ltime =  systime(1)
mtime =  ltime
chi2 = -1
print,  'Press enter to end...'
while (1) do begin
 while (mtime-ltime le  0) do begin
  openr,  1,  'thread00.fits',  error = err
  if (err eq 0) then begin
   stat =  fstat(1)
   close,  1
   mtime =  stat.mtime
 endif else print,  'Trouble opening file (ignore if occasional message only)'
  wait,  0.3
  if (get_kbrd(0) ne '') then goto,  finish
 endwhile
 ltime = mtime
 wset,  0
 im =  readfits('thread00.fits',  head)
 image_cont,  im^0.5,  /noc
 wset,  1
 newchi2 = sxpar(head,  'CHISQR')
 if (chi2[0] eq -1) then chi2 = newchi2 else begin
  chi2 =  [chi2, newchi2]
  plot,  chi2,  yr = [0, min([min(chi2)*2,  max(chi2)])]
 endelse
endwhile
finish:
print,  'Enter pressed, exiting...'

end
