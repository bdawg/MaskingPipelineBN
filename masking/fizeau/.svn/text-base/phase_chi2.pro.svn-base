;A function for phase chi^2 in the presence of wrapping
function phase_chi2, p
common phasestuff, fitmat, ph_mn, ph_err
return, total(modsq(1.-exp(complex(0,[ph_mn,0] - reform(p#fitmat))))/[ph_err,0.01]^2)
;return, total(modsq(1.-exp(complex(0,ph_mn - reform(p#fitmat))))/ph_err^2)
end
