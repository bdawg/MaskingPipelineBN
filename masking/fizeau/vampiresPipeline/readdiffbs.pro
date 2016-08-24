; Reads sets of bs___.idlvar files and extracts v2s and bispectra.
; Called by diffcal_vampires.script.
;
; Arrays are of form [bls, frames, Wollaston, LCVR]

pro readdiffbs, prefix, hwp0FileNums, hwp225FileNums, hwp45FileNums, hwp675FileNums, $
                h0_v2s_all, h0_bss_all, h225_v2s_all, h225_bss_all, $
                h45_v2s_all, h45_bss_all, h675_v2s_all, h675_bss_all, $
                u_coords, v_coords


; Zero v2s variables
h0_v2s=0
h225_v2s=0
h45_v2s=0
h675_v2s=0


; HWP=0 set:
nfiles=n_elements(hwp0FileNums)
for readFiles = 0,nfiles-1 do begin    
    for wollsChan = 0,1 do begin
        for lcvrState = 0,1 do begin
            if (wollsChan eq 0) AND (lcvrState eq 0) then fileSuf = '_1_A'
            if (wollsChan eq 0) AND (lcvrState eq 1) then fileSuf = '_1_B'
            if (wollsChan eq 1) AND (lcvrState eq 0) then fileSuf = '_2_A'
            if (wollsChan eq 1) AND (lcvrState eq 1) then fileSuf = '_2_B'

            filename='bs_'+prefix+strn(hwp0FileNums[readFiles],len=4,padchar='0')+ $
              fileSuf+'.idlvar'
            restore,filename
            if n_elements(h0_v2s) eq 1 then begin
                nfrms = (size(v2_all))[1]
                nbls = (size(v2_all))[2]
                ncps = (size(cp))[1]
                h0_v2s=dblarr(nbls,nfrms,2,2)
                h0_bss=dcomplexarr(ncps,nfrms,2,2)
                h0_u=dblarr(nbls,2,2)
                h0_v=dblarr(nbls,2,2)
            endif
            h0_v2s[*,*,wollsChan,lcvrState]=transpose(v2_all)
            h0_bss[*,*,wollsChan,lcvrState]=transpose(bs_all)
            h0_u[*,wollsChan,lcvrState]=u
            h0_v[*,wollsChan,lcvrState]=v
        endfor
    endfor

    if readFiles eq 0 then begin
        h0_v2s_all=h0_v2s
        h0_bss_all=h0_bss
        h0_u_all=h0_u
        h0_v_all=h0_v
        h0_v2s=0
    endif else begin
        h0_v2s_all=[[h0_v2s_all],[h0_v2s]]
        h0_bss_all=[[h0_bss_all],[h0_bss]]
        h0_v2s=0
    endelse

 endfor

; HWP=22.5 set:
nfiles=n_elements(hwp225FileNums)
for readFiles = 0,nfiles-1 do begin    
    for wollsChan = 0,1 do begin
        for lcvrState = 0,1 do begin
            if (wollsChan eq 0) AND (lcvrState eq 0) then fileSuf = '_1_A'
            if (wollsChan eq 0) AND (lcvrState eq 1) then fileSuf = '_1_B'
            if (wollsChan eq 1) AND (lcvrState eq 0) then fileSuf = '_2_A'
            if (wollsChan eq 1) AND (lcvrState eq 1) then fileSuf = '_2_B'

            filename='bs_'+prefix+strn(hwp225FileNums[readFiles],len=4,padchar='0')+ $
              fileSuf+'.idlvar'
            restore,filename
            if n_elements(h225_v2s) eq 1 then begin
                nfrms = (size(v2_all))[1]
                nbls = (size(v2_all))[2]
                ncps = (size(cp))[1]
                h225_v2s=dblarr(nbls,nfrms,2,2)
                h225_bss=dcomplexarr(ncps,nfrms,2,2)
                h225_u=dblarr(nbls,2,2)
                h225_v=dblarr(nbls,2,2)
            endif
            h225_v2s[*,*,wollsChan,lcvrState]=transpose(v2_all)
            h225_bss[*,*,wollsChan,lcvrState]=transpose(bs_all)
            h225_u[*,wollsChan,lcvrState]=u
            h225_v[*,wollsChan,lcvrState]=v
        endfor
    endfor

    if readFiles eq 0 then begin
        h225_v2s_all=h225_v2s
        h225_bss_all=h225_bss
        h225_u_all=h225_u
        h225_v_all=h225_v
        h225_v2s=0
    endif else begin
        h225_v2s_all=[[h225_v2s_all],[h225_v2s]]
        h225_bss_all=[[h225_bss_all],[h225_bss]]
        h225_v2s=0
    endelse

endfor

; HWP=45 set:
nfiles=n_elements(hwp45FileNums)
for readFiles = 0,nfiles-1 do begin    
    for wollsChan = 0,1 do begin
        for lcvrState = 0,1 do begin
            if (wollsChan eq 0) AND (lcvrState eq 0) then fileSuf = '_1_A'
            if (wollsChan eq 0) AND (lcvrState eq 1) then fileSuf = '_1_B'
            if (wollsChan eq 1) AND (lcvrState eq 0) then fileSuf = '_2_A'
            if (wollsChan eq 1) AND (lcvrState eq 1) then fileSuf = '_2_B'

            filename='bs_'+prefix+strn(hwp45FileNums[readFiles],len=4,padchar='0')+ $
              fileSuf+'.idlvar'
            restore,filename
            if n_elements(h45_v2s) eq 1 then begin
                nfrms = (size(v2_all))[1]
                nbls = (size(v2_all))[2]
                ncps = (size(cp))[1]
                h45_v2s=dblarr(nbls,nfrms,2,2)
                h45_bss=dcomplexarr(ncps,nfrms,2,2)
                h45_u=dblarr(nbls,2,2)
                h45_v=dblarr(nbls,2,2)
            endif
            h45_v2s[*,*,wollsChan,lcvrState]=transpose(v2_all)
            h45_bss[*,*,wollsChan,lcvrState]=transpose(bs_all)
            h45_u[*,wollsChan,lcvrState]=u
            h45_v[*,wollsChan,lcvrState]=v
        endfor
    endfor

    if readFiles eq 0 then begin
        h45_v2s_all=h45_v2s
        h45_bss_all=h45_bss
        h45_u_all=h45_u
        h45_v_all=h45_v
        h45_v2s=0
    endif else begin
        h45_v2s_all=[[h45_v2s_all],[h45_v2s]]
        h45_bss_all=[[h45_bss_all],[h45_bss]]
        h45_v2s=0
    endelse

endfor

; HWP=67.5 set:
nfiles=n_elements(hwp675FileNums)
for readFiles = 0,nfiles-1 do begin    
    for wollsChan = 0,1 do begin
        for lcvrState = 0,1 do begin
            if (wollsChan eq 0) AND (lcvrState eq 0) then fileSuf = '_1_A'
            if (wollsChan eq 0) AND (lcvrState eq 1) then fileSuf = '_1_B'
            if (wollsChan eq 1) AND (lcvrState eq 0) then fileSuf = '_2_A'
            if (wollsChan eq 1) AND (lcvrState eq 1) then fileSuf = '_2_B'

            filename='bs_'+prefix+strn(hwp675FileNums[readFiles],len=4,padchar='0')+ $
              fileSuf+'.idlvar'
            restore,filename
            if n_elements(h675_v2s) eq 1 then begin
                nfrms = (size(v2_all))[1]
                nbls = (size(v2_all))[2]
                ncps = (size(cp))[1]
                h675_v2s=dblarr(nbls,nfrms,2,2)
                h675_bss=dcomplexarr(ncps,nfrms,2,2)
                h675_u=dblarr(nbls,2,2)
                h675_v=dblarr(nbls,2,2)
            endif
            h675_v2s[*,*,wollsChan,lcvrState]=transpose(v2_all)
            h675_bss[*,*,wollsChan,lcvrState]=transpose(bs_all)
            h675_u[*,wollsChan,lcvrState]=u
            h675_v[*,wollsChan,lcvrState]=v
        endfor
    endfor

    if readFiles eq 0 then begin
        h675_v2s_all=h675_v2s
        h675_bss_all=h675_bss
        h675_u_all=h675_u
        h675_v_all=h675_v
        h675_v2s=0
    endif else begin
        h675_v2s_all=[[h675_v2s_all],[h675_v2s]]
        h675_bss_all=[[h675_bss_all],[h675_bss]]
        h675_v2s=0
    endelse

endfor

u_coords = h0_u[*,0,0]
v_coords = h0_v[*,0,0]

end
