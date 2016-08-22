function vartype, variable, sz, CODE=code
;+
; NAME:
;	vartype
;
; PURPOSE:
;	Check the type of IDL variables.
;
; CALLING:
;	type_string = vartype( variable )
; or:
;	type_code = vartype( variable, /CODE )
;
; INPUTS:
;	variable = anything.
;
; KEYWORDS:
;	/CODE : causes the integer IDL type code to be returned,
;		instead of a string describing variable type.
; OUTPUTS:
;	sz = (optional) the result of the function size( variable ).
;
;	Function returns string describing variable type,
;	or if /CODE then the integer IDL type code (0 to 8).
;
; HISTORY:
;	Written, Frank Varosi NASA/GSFC 1989.
;-
	sz = size( variable )
	type = sz( sz(0)+1 )
	if keyword_set( code ) then return,type

	CASE type OF
	1:	typename = "BYTE"
	2:	typename = "INTEGER SHORT"
	3:	typename = "INTEGER LONG"
	4:	typename = "FLOATING"
	5:	typename = "FLOATING DOUBLE"
	6:	typename = "FLOATING COMPLEX"
	7:	typename = "STRING"
	8:	typename = "STRUCTURE"
	else:	typename = ""
	ENDCASE

return, typename
end
