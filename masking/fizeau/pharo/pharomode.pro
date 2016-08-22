function pharomode, head

; take a header and make a modestring
;
;
  
    t_int = sxpar(head,'T_INT')
    writemod =sxpar(head,'WRITEMOD')
    nendptfr = sxpar(head,'NENDPTFR')
    npausefr = sxpar(head,'NPAUSEFR')
    naxis1 = sxpar(head,'NAXIS1')
    naxis2 = sxpar(head,'NAXIS2')
    
    mode = strc(t_int)+':' + $
           strc(nendptfr) +'+'+ strc(npausefr) +'+'+ strc(nendptfr) + ':' +$
           strc(naxis1) + 'x' + strc(naxis2) +':' + $
           strc(writemod)
    return, mode
end
