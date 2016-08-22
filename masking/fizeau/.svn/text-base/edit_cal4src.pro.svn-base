;######## RELEASE CLIP HERE ######## 
;
; edit_cal4src.pro
; Sep 2009
; Paul Stewart
; 
; Takes input of a cal4src array and an array of lables (filenumbers perhaps)
; Allows manual alterations to be made to the matrix
; Outputs the ammended matrix
;

function edit_cal4src, cal4src, list
if n_params() ne 2 then begin
    print, " ### ERROR ### Invalid Input Params to cal4src editor"
    return, 0
endif

xlist = list
ylist = reverse(xlist)

if max(cal4src) gt 1 or min(cal4src) lt 0 then begin
  print, " ### ERROR ### Invalid Input array to cal4src editor"
  return, 0
endif

check = 0
;read, "For whole matrix enter 0, for active quadrant only enter 1: ", check
if check eq 1 then begin
    master = cal4src
    w = where(cal4src eq 1)
    cal4src = cal4src[w]
    cal4src = reform(cal4src,fix(sqrt(n_elements(w))),fix(sqrt(n_elements(w))))
    
endif

;set up display array
cal4src = reverse(cal4src,2)
s = size(cal4src)
sx = s[1]-1
sy = s[2]-1
sc = 480 ;scale number to fit into 480 pixels
gridx = indgen(s[1]+1) * (sc/s[1])
gridy = indgen(s[2]+1) * (sc/s[2])
textx = (findgen(s[1])/s[1]*1.)*.8+.1+(.4/s[1]*1.)
texty = (findgen(s[2])/s[2]*1.)*.8+.1+(.4/s[2]*1.)

;Setting up the display
loadct, 39
WINDOW,0,XSIZE=600,YSIZE=600, title="cal4src Matrix Editor"
x = 0

;Main Loop
while x ge 0 do begin
    ;generate output and display
    display = 150*congrid(cal4src,sc,sc)
    display[gridx,*]=255
    display[*,gridy]=255
    tv, display, .1, .1, /normal
    xyouts,.5,.97,"cal4src Matrix Editor",/normal, align=.5
    xyouts,.5,.94,"Click outside the grid to accept and exit",/normal, align=.5
    xyouts, textx, .07, xlist,/normal,align=.5
    xyouts, .09, texty, ylist,/normal,align=1

    cursor, x, y, /normal,/down
    if x lt 0.1 or y lt .1 or x gt .9 or y gt .9 then break
    x = floor(((x-.1)/.8)*(s[1]))
    y = floor(((y-.1)/.8)*(s[2]))
    if cal4src[x,y] eq 0 then cal4src[x,y] = 1 else cal4src[x,y] = 0
endwhile
wdelete

cal4src = reverse(cal4src,2)

;reinsert array into master, then rename to cal4src
if check eq 1 then begin
    master(w) = cal4src
    cal4src = master
endif


print, cal4src
return, cal4src
end

