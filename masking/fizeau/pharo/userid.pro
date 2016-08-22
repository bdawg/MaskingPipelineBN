function userid

;+
; return string with user and host, padded by spaces
; 03/15/99 MCL
;
; make sure to return a scalar, not a 1-element array (stupid IDL)
; 10/05/00 MCL
;-


spawn, 'echo $USER@$HOST', id

return, ' '+id(0)+' '

end

