
;; This removes the spaces and dashes in the target name and replaces with
;; underscores

FUNCTION FORMAT_TARGET_NAME, name

   ; Convert Target Name spaces to underscores
   TargetByte = byte( Name )
   SpaceByte = (byte( " " ))[0]
   DashByte  = (byte( "-" ))[0]
   PeriodByte= (byte( "." ))[0]
   PlusByte  = (byte( "+" ))[0]
   
   UnderByte = (byte( "_" ))[0]
   w1 = Where( TargetByte EQ SpaceByte, C1 )
   w2 = Where( TargetByte EQ DashByte, C2 )
   w3 = Where( TargetByte EQ PeriodByte, C3 )
   w4 = Where( TargetByte EQ PlusByte, C4 )
   
   ; Do the swap: " " for "_"
   If( C1 NE 0 ) THEN TargetByte[ w1 ] = UnderByte
   If( C2 NE 0 ) THEN TargetByte[ w2 ] = UnderByte
   If( C3 NE 0 ) THEN TargetByte[ w3 ] = UnderByte
   If( C4 NE 0 ) THEN TargetByte[ w4 ] = UnderByte

   NewName = String( TargetByte )

   Return, NewName
END
