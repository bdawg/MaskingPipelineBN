;;ALL FEAR THE MIGHTY MEGA CAL4SRC!
;;
;;Designed to be a single program to replace all options for
;;cal4src. Recognises blocks from the tsize variable.
;;
;;Allows the user edit cal4src matrices in the following ways:
;;
;;1. Select which observing blocks should be used by clicking on them.
;;
;;2. Select individual files that are part of each block
;;
;;3. Use match quads (Only calibrate files if the star is in the same
;;quadrant of the chip.
;;
;;4. Select All
;;
;;5. Select None
;;
;;9. Allow selection of calibrators-calibrator blocks
;;
;;Note that when the cal4src array is displayed, [0,0] will appear
;;white, so that the rest of the array can appear green
;;
;;Set /unlock_all to calibrate cals with other cals. Will allow any
;;block to be calibrated with any other block
;;
function mega_cal4src,cubeinfo,unlock_all=unlock_all
old_pmulti=!p.multi
!p.multi=0

;;first, get the old cal4src
restore,cubeinfo
cal4src=olog.cal4src
tsize=olog.cube_tsize
new=0*cal4src
nfiles=n_elements(cal4src[*,0])

;;find which files are calibrators and which are sources
cals=where(tsize ge 0,complement=src)
nc=n_elements(cals)
ns=n_elements(src)

;;find where the files change from cals to sources. EDIT: Now uses
;;changes in tsize to find block changes. So two targets observed in a
;;row show up as separate blocks.
;breaks=where(cals[1:*]-cals[*] ne 1)
;breaks2=where(src[1:*]-src[*] ne 1)
if n_elements(cals) ne 1 then begin
   breaks=where(tsize[cals[1:*]]-tsize[cals[*]] ne 0)
   ncals=n_elements(breaks)+1
endif else ncals=1
if n_elements(src) ne 1 then begin
   breaks2=where(tsize[src[1:*]]-tsize[src[*]] ne 0)
   nsrc=n_elements(breaks2)+1
endif else nsrc=1

n_blocks=nsrc+ncals

;;make an array with the start and end of each calibrator and block.
blocks=fltarr(2,ncals)
ix1=cals[0]
if ncals eq 1 then begin
   ix2=cals
   ix1=cals
   blocks[0]=ix1
   blocks[1]=ix2
endif else begin
   for i =0,ncals-2 do begin
      ix2=cals[breaks[i]]
      blocks[0,i]=ix1
      blocks[1,i]=ix2
      ix1=cals[breaks[i]+1]    
   endfor
   ;;but this will miss the last block, so include it manually
   blocks[0,i]=ix1
   blocks[1,i]=cals[nc-1]
endelse

;;now do the same thing for the sources
blocks2=fltarr(2,nsrc)
ix1=src[0]
if nsrc eq 1 then begin
   ix2=src
   ix1=src
   blocks2[0]=ix1
   blocks2[1]=ix2
endif else begin
   for i =0,nsrc-2 do begin
      ix2=src[breaks2[i]]
      blocks2[0,i]=ix1
      blocks2[1,i]=ix2
      ix1=src[breaks2[i]+1]    
   endfor
   blocks2[0,i]=ix1
   blocks2[1,i]=src[ns-1]
endelse

if (nsrc +ncals) lt 10 then begin
   print,"It looks like you don't have many observations. It will be easier to use option 2!"
endif

;;Label the cals and the targets:
;;Get ra and dec
ra=olog.ra
dec=olog.dec

cal_radec=fltarr(2,ncals)
src_radec=fltarr(2,nsrc)
cal_label=strarr(ncals)
;src_label=replicate("Target",nsrc);assume there is one target
src_label=strarr(nsrc)

;;fill in radec arrays
for i=0,ncals-1 do begin
    cal_radec[0,i]=ra[blocks[i]]
    cal_radec[1,i]=dec[blocks[i]]
endfor
for i=0,nsrc-1 do begin
    src_radec[0,i]=ra[blocks2[i]]
    src_radec[1,i]=dec[blocks2[i]]
endfor

;;name the stars:
;;check none are targets!
for i=0,nsrc-1 do begin
   w=where(sqrt(total((rebin(src_radec[*,i],[2,ncals])-cal_radec)^2,1)) lt 0.01)
   cal_label[w]='T'+strcompress(string(i),/remove_all)+'?'
endfor
;;Cals
tar=0
for i=0,ncals-1 do begin
   if not keyword_set(cal_label[i]) then begin
      tar+=1
      cal_label[i]='C'+strcompress(string(tar),/remove_all)
      ;;now assume the rest of the blocks within 10 arcsec are of
      ;;the same cal
      w=where(sqrt(total((rebin(cal_radec[*,i],[2,ncals])-cal_radec)^2,1)) lt 0.01)
      cal_label[w]='C'+strcompress(string(tar),/remove_all)
   endif
endfor

;;Targets:
tar=0
for  i=0,nsrc-1 do begin
   if not keyword_set(src_label[i]) then begin
      tar+=1
      src_label[i]='T'+strcompress(string(tar),/remove_all)
      ;;now assume the rest of the blocks within 10 arcsec are of
      ;;the same target
      w=where(sqrt(total((rebin(src_radec[*,i],[2,ncals])-src_radec)^2,1)) lt 0.01)
      src_label[w]='T'+strcompress(string(tar),/remove_all)
   endif
endfor

;;set up a few numbers for plotting. These are based on what looks good on ACC's screen
yoffset=-nfiles*0.05;;for src?
xoffset=-nfiles*0.1;;for cals


;;what are the x and y values needed so that everything sits between
sz=size(cal4src)
xvals=sz[1]*findgen(sz[1])/(sz[1]-1)+0.5
yvals=sz[1]*findgen(sz[2])/(sz[1]-1)+0.5

select_option: ;;we will loop from here down, so multiple options can be selected

;;We want to ask the user what they want to do, so prompt them.
print,'Select an option'
read,option,prompt='(0) Accept. (1) Select blocks. (2) Select files. (3) Match quads. (4) Select All. (5) Select None (9) Allow Cal-Cal blocks:'
if option eq 1 then begin
    print,'Select which block you want. Possible blocks are white, selected blocks are green'
    print,'Click outside the plot when finished'
    ;;Allow the user to select as many blocks the want. Clicking outside
    ;;the plot area will exit the loop.
    while (1) do begin
        ;;ask which block they want:
        plot_im=2.*cal4src-new
        plot_im[0,0]=2 ;;so that everything else is green.
        image_cont,plot_im,/n,tit='Current cal4src',xval=xvals,yval=yvals
        xyouts,blocks[0,*],replicate(yoffset,ncals),cal_label,/data,charsize=2,charthick=2
        xyouts,replicate(xoffset,nsrc),((blocks2[0,*]+blocks2[1,*])/2.),src_label,/data,charsize=2,charthick=2
        cursor,x,y & wait,0.2
        
        ;;round down, since indices are integers while x and y are not:
        x=floor(x-0.5)
        y=floor(y-0.5)
        
        if x lt 0 or y lt 0 then goto, select_option ;;exits the loop, going back to the select option stage
        
        ;;what block did you select?
        wx=where(blocks[1,*] ge x and blocks[0,*] le x)
        wy=where(blocks2[1,*] ge y and blocks2[0,*] le y)
        
        if n_elements(wx)+n_elements(wy) ne 2 then stop ;;this should have flagged one block in x and one block in y
        if wx[0] eq -1 or wy[0] eq -1 then begin
            print,'You missed! Try again!'
            goto,try_again_blocks
        endif;;for some reason it did not find the block
        
        cix1=blocks[0,wx[0]]
        cix2=blocks[1,wx[0]]
        six1=blocks2[0,wy[0]]
        six2=blocks2[1,wy[0]]
        
        
        if median(new[cix1:cix2,six1:six2]) eq 0 then new[cix1:cix2,six1:six2]=1 else if median(new[cix1:cix2,six1:six2]) eq 1 then new[cix1:cix2,six1:six2]=0
        
        try_again_blocks:
        
        ;;The unlock_all keyword allows cal-cal blocks to be
        ;;selected. This allows calibrators to be considered as targets:
        if (keyword_set(unlock_all) and (wx[0] eq -1 or wy[0] eq -1) )then begin
            if wx[0] eq -1 then wx=where(blocks[1,*] ge x and blocks[0,*] le x)
            if wy[0] eq -1 then wy=where(blocks[1,*] ge y and blocks[0,*] le y)
            ;;take both indices from the calibrator blocks
            cix1=blocks[0,wx[0]]
            cix2=blocks[1,wx[0]]
            six1=blocks[0,wy[0]]
            six2=blocks[1,wy[0]]
            if median(new[cix1:cix2,six1:six2]) eq 0 then new[cix1:cix2,six1:six2]=1 else if median(new[cix1:cix2,six1:six2]) eq 1 then new[cix1:cix2,six1:six2]=0
        endif
        wait,0.1
    endwhile
endif else if option eq 2 then begin
    ;;now, allow the user to select files individually:
    print,'Click to select/deselect individual files. Possible files are white, selected files are green.'
    print,'Click outside the plot when finished'
    ;;Allow the user to select as many blocks the want. Clicking outside
    ;;the plot area will exit the loop.
    while (1) do begin
        ;;ask which block they want:
        plot_im=2.*cal4src-new
        plot_im[0,0]=2 ;;so that everything else is green.
        image_cont,plot_im,/n,tit='Current cal4src',xval=xvals,yval=yvals
        xyouts,blocks[0,*],replicate(yoffset,ncals),cal_label,/data,charsize=2,charthick=2
        xyouts,replicate(xoffset,nsrc),((blocks2[0,*]+blocks2[1,*])/2.),src_label,/data,charsize=2,charthick=2
        cursor,x,y & wait,0.2

        ;;round down, since indices are integers while x and y are not:
        x=floor(x-0.5)
        y=floor(y-0.5)

        if x lt 0 or y lt 0 then goto, select_option   

        if x eq y then begin
           print,'Warning! You selected a block to calibrate with itself!'
        endif ;;This might be useful to know...
        
        if new[x,y] eq 1 then new[x,y]=0 else if new[x,y] eq 0 then new[x,y]=1
        
        try_again_indiv:
        wait,0.1
    endwhile
endif else if option eq 3 then begin
    ;;matchquads: only calibrate files if they are from the same quadrant
    quadmatch = intarr(nfiles,nfiles)
    for xix = 0, nfiles-1 do for yix = 0, nfiles-1 do if olog.quad[xix] eq olog.quad[yix] then quadmatch[xix,yix]=1
    new*=quadmatch
endif else if option eq 4 then begin
    ;;Select all!
    new*=0
    for i=0,n_elements(src)-1 do new[cals,src[i]]=1
endif else if option eq 5 then begin
    ;;select none!
    new*=0
endif else if option eq 0 then begin
    ;;Done! Exit loop
    goto, done
 endif else if option eq 9 then begin
    if keyword_set(unlock_all) then unlock_all=0 else unlock_all=1
    if keyword_set(unlock_all) then print,'unlocked calibrator-calibrator blocks' else print,'locked calibrator-calibrator blocks'
 endif else if option eq stop then begin
    stop
 endif else begin
    print,'Incorrect option. Try again!'
 endelse   
;;display the resulting cal4src:
plot_im=2.*cal4src-new
plot_im[0,0]=2 ;;so that everything else is green.
image_cont,plot_im,/n,tit='Current cal4src',xval=xvals,yval=yvals
xyouts,blocks[0,*],replicate(yoffset,ncals),cal_label,/data,charsize=2,charthick=2
xyouts,replicate(xoffset,nsrc),((blocks2[0,*]+blocks2[1,*])/2.),src_label,/data,charsize=2,charthick=2

;;now send back to the start, so we can continue editing.
goto,select_option

done:
image_cont,new,/n,tit='Final cal4src Matrix',xval=xvals,yval=yvals
;xyouts,blocks[0,*],replicate(yoffset,ncals),cal_label,/data,charsize=2,charthick=2
;xyouts,replicate(xoffset,nsrc),((blocks2[0,*]+blocks2[1,*])/2.),src_label,/data,charsize=2,charthick=2

!p.multi=old_pmulti

return,new
end
