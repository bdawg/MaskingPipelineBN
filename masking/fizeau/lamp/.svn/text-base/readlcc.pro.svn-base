;This function will read the LUCKY format compressed fits files.
function readlcc,  filename,  head,  silent = silent
 
 ;Uncompress the file. !!! Needs lucky2fits in the PATH
 spawn, 'lucky2fits '+filename, unit=unit
 a = readfits(unit,  head,  silent = silent)
 free_lun,  unit
 return,  a

end
