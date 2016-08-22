;;Edits cal4src, so that only calibrators observed close in time to
;;each source will be used for calibration.
;;
;;Inputs:
;;        radius: How many cals either side of each source to use.
;;        cubeinfo: the filename of the cubeinfo file
;;        
;;Options:
;;        all: Use all cals with all sources
;;
function cal4src_diag,cubeinfo,radius,all=all
restore,cubeinfo
cal4src=olog.cal4src
new=0*cal4src

nfiles=n_elements(cal4src[*,0])
if keyword_set(all) then radius= nfiles+1 ;;a number that is bigger than the number of blocks.

;;find which files are calibrators and which are sources
cals=where(olog.cube_tsize ge 0,complement=src)
nc=n_elements(cals)
ns=n_elements(src)

;;find where the files change from cals to sources
breaks=where(cals[1:*]-cals[*] ne 1)
breaks2=where(src[1:*]-src[*] ne 1)
ncals=n_elements(breaks)+1
nsrc=n_elements(breaks2)+1

;;make an array with the start and end of each calibrator and block.
blocks=fltarr(2,ncals)
ix1=cals[0]
for i =0,ncals-2 do begin
    ix2=cals[breaks[i]]
    blocks[0,i]=ix1
    blocks[1,i]=ix2
    ix1=cals[breaks[i]+1]    
endfor
;;but this will miss the last block, so include it manually
blocks[0,i]=ix1
blocks[1,i]=cals[nc-1]


blocks2=fltarr(2,nsrc)
ix1=src[0]
for i =0,nsrc-2 do begin
    ix2=src[breaks2[i]]
    blocks2[0,i]=ix1
    blocks2[1,i]=ix2
    ix1=src[breaks2[i]+1]    
endfor
blocks2[0,i]=ix1
blocks2[1,i]=src[ns-1]

;;now go through each cal, and make cal4src for it
for i=0,nsrc-1 do begin

    six1=blocks2[0,i]
    six2=blocks2[1,i]

    for j=0,ncals-1 do begin

        cix1=blocks[0,j]
        cix2=blocks[1,j]
        
        ;;this will be different depending on whether a source or
        ;;cal is observed first:
        if cals[0] eq 0 then begin;;if a cal is first
            if ((j-i) le radius) and ((j-i) ge -(radius-1)) then  new[cix1:cix2,six1:six2]=1
        endif else if src[0] eq 0 then begin
            if ((j-i) le (radius-1)) and ((j-i) ge -radius) then new[cix1:cix2,six1:six2]=1
        endif else stop ;;something went wrong... (neither cal nor source is first?)

        
    endfor
endfor

;;update the new cal4src, and save the updated cubeinfo file
;clog.cal4src=new
;save,olog,plog,clog,framestats,quan,calib_src_quan,u_save,v_save,u1,v1,src,filename=cubeinfo

return,new
end
