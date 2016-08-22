
;; This function sorts an array of structures by the values stored in
;; the tags listed in "order"  (see example below)  The structures are
;; ordered first by the first element in the ORDER variable and ties
;; are broken by ordering the second element in the ORDER variable,
;; and so on
;;
;; Struct = replicate( { A: 0L, B: 0L }, 10 )
;; Order  = [ "B", "A" ] or [ "A", "B" ]
;;
;;
;;


function Sort_Array_Of_Structs_By_Tags, Struct, Order, ReturnStruct_in

  If( Keyword_Set(ReturnStruct_in) EQ 0 ) THEN ReturnStruct_in = 0
  ReturnStruct = ReturnStruct_in

  For TagIndex = n_elements( Order ) - 1, 0, -1 DO BEGIN

     StructIndex = Where( Tag_Names( Struct ) EQ StrUpCase( Order[TagIndex] ) )
     Values = Struct.( StructIndex[0] )

     s = BSort( Values )

     NewStruct = Struct
     For i = 0, n_elements( s ) - 1 DO NewStruct[i] = Struct[ s[i] ]
     Struct = NewStruct

  ENDFOR

  If( ReturnStruct ) THEN Return, Struct ELSE Return, s

END
