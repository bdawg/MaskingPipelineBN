
;; -----------------------------------------------------------------------------
;;                           FIT A BINARY SYSTEM
;; -----------------------------------------------------------------------------

pro analyse_binary, datafile, parmlimits, return_info, use_dwb_errors = use_dwb_errors

;; --------------------------------------------------------
;; Extract closures from the file
;; 
;;  If using DWB's method, grab data from correct file.
;; --------------------------------------------------------

  extract_t3data, file = datafile, t3data

  print, 'Median closure phase error: ', median(t3data.t3phierr),' degrees' 

  t3null = t3data
  t3null.t3phi = t3data.t3phi * 0.0

;; --------------------------
;; Set Parm Limits
;; --------------------------

  minSep = ParmLimits[0].min
  maxSep = ParmLimits[0].max
  nSep   = ParmLimits[0].n

  minAz = ParmLimits[1].min
  maxAz = ParmLimits[1].max
  nAz   = ParmLimits[1].n

  minCRat = ParmLimits[2].min
  maxCRat = ParmLimits[2].max
  nCRat   = ParmLimits[2].n

  Sep  = minSep + (maxSep - minSep) * findgen(nSep)/nSep
  Az   = minAz  + (maxAz - minAz) * findgen(nAz) / nAz
  Crat = 10.^( minCrat + (maxCrat - minCrat) * findgen(nCrat)/nCrat )

  chi2cube = fltarr(nSep, nAz, nCrat)
  nullchi2 = total( (t3data.t3phi/t3data.t3phierr)^2 )
  for i = 0, nSep-1 do begin
     for j = 0, nAz - 1 do begin
        for k = 0, nCrat - 1 do begin
           p0 = [Sep[i], Az[j], Crat[k]]
           
           ;; ---- compute vis & closures ----
           t3model   = binary_t3data_fast(p0, t3data=t3data)
           
           ;; ---- evaluate the chi2 ----
           chi2 = binary_chi2_fast(p0, t3model=t3model, t3data=t3data)
           chi2cube[i,j,k] = chi2
        Endfor
     Endfor
     print, i
  Endfor


;; ----------------------------------------------
;;             Display the chi2 in the  
;;               Sep-Az plane
;; ----------------------------------------------

;  Window, 2, xsize=750, ysize=750
  loadct, 0

  CratNum = 0
  PlotType = 1

  minall = fltarr(nSep,nAz)
  for i = 0, nSep-1 do begin
     for j = 0, nAz - 1 do begin
        minall[i,j] = MIN( chi2cube[i,j,*] )
     ENDFOR
  ENDFOR

  REPEAT BEGIN
     If (PlotType EQ 1 ) THEN BEGIN
        xyarray = reform(chi2cube[*,*, CratNum])
        title = 'Chi2 (Con. Rat = ' + String( Crat[CratNum] ) + ')'
     ENDIF  

     If (PlotType EQ 2 ) THEN BEGIN
        xyarray = minall
        title = 'Chi2 (Con. Rat = ' + String( Crat[CratNum] ) + ')'
     ENDIF  
     
     If (PlotType EQ 3 ) THEN BEGIN
        xyarray = minall < nullchi2
        title = 'Min Chi2 (Crat), but nothing above NullChi2'
     ENDIF

     contour, xyarray, Sep, Az, /fill, tit=title , nlevel=20, xtitle='angular sep (mas)', ytitle = 'azimut (degrees)', xrange=[minSep, maxSep], yrange=[minAz, maxAz]
     contour, xyarray, Sep, Az, tit='chi2', nlevel=20, /overplot, /follow, xrange=[minSep, maxSep], yrange=[minAz, maxAz]
     
     Key = GET_KBRD(/KEY_NAME)

     If ( Key EQ '1' ) THEN PlotType = 1
     If ( Key EQ '2' ) THEN PlotType = 2
     If ( Key EQ '3' ) THEN PlotType = 3
     If ( Key EQ '.' ) THEN CratNum = MIN([CratNum + 1, nCrat-1]) 
     If ( Key EQ ',' ) THEN CratNum = MAX([CratNum - 1, 0 ])
     If ( Key EQ 'q' ) THEN BEGIN
        Print, 'Click on Maximum'
        cursor, optiSep, optiAz
        
        Print, 'Optimum Separation: ', optiSep, ' mas at ', optiAz, ' degrees.'
        BREAK
     ENDIF
  ENDREP UNTIL (1 EQ 0)  ;; Forever

;; --------------------------------------------------------
;; Display the chi2 as a function of Contrast Ratio
;; ---------------------------------------------------------

; Increase Crat resolution
  nCrat = nCrat*3

 Crat = 10.^( minCrat + (maxCrat - minCrat) * findgen(nCrat)/nCrat )

  chi2vCrat = fltarr(nCrat)
  For i = 0, nCrat - 1 do BEGIN
     p0 = [optiSep, optiAz, Crat[i]]
     
     ;; ---- compute vis & closures ----
     t3model   = binary_t3data_fast(p0, t3data=t3data)
     
     ;; ---- evaluate the chi2 ----
     chi2 = binary_chi2_fast(p0, t3model=t3model, t3data=t3data)
     chi2vCrat[i] = chi2
  Endfor

  maxyrange = Min(chi2vCrat) * 2.5
  plot, Crat, chi2vCrat, xtitle = "Contrast Rato", ytitle = "Chi2", yrange = [0, maxyrange], /xlog

  REPEAT BEGIN
     Key = GET_KBRD(/KEY_NAME)
  ENDREP UNTIL (Key EQ 'q')


  Print, 'Click on Minimum'
  cursor, optiCrat, y1

  print, 'Optimum Contrast Ratio = ', optiCrat


;; -----------------------------------------------------------------
;; Zoom in on Sep-Az plane at chosen Crat
;; -----------------------------------------------------------------

  CratPRange = .30
  
  SepRange = ( MaxSep - MinSep ) / 3.
  AzRange  = ( MaxAz - MinAz ) / 3.

  ZoomSep  = OptiSep + SepRange * ( findgen(nSep)/nSep - .5 )
  ZoomAz   = OptiAz  + AzRange * ( findgen( nAz ) / nAz - .5 )
  ZoomCrat = optiCrat

  chi2cube = fltarr(nSep, nAz)

  for i = 0, nSep-1 do begin
     for j = 0, nAz - 1 do begin
        p0 = [ZoomSep[i], ZoomAz[j], optiCrat]
        
        ;; ---- compute vis & closures ----
        t3model   = binary_t3data_fast(p0, t3data=t3data)
        
        ;; ---- evaluate the chi2 ----
        chi2 = binary_chi2_fast(p0, t3model=t3model, t3data=t3data)
        
        chi2cube[i,j] = chi2    
     Endfor
     print, i
  Endfor

;  Window, 2, xsize=750, ysize=750
  !p.multi=0
  loadct, 0
  contour, Chi2cube, ZoomSep, ZoomAz, /fill, tit='chi2', nlevel=20, xtitle='angular sep (mas)', ytitle = 'azimut (degrees)', xrange=[minSep, maxSep], yrange=[minAz, maxAz]
  contour, Chi2Cube, ZoomSep, ZoomAz, tit='chi2', nlevel=20, /overplot, /follow, xrange=[minSep, maxSep], yrange=[minAz, maxAz]

  REPEAT BEGIN
     Key = GET_KBRD(/KEY_NAME)
  ENDREP UNTIL (Key EQ 'q')

  Print, 'Click on Maximum'
  cursor, optiSep, optiAz

  Print, 'Optimum Separation: ', optiSep, ' mas at ', optiAz, ' degrees at CRatio ', optiCrat, '.'


;; ------------------------------------------------------------
;;  Use mpfit and the OptiParameters as starting points to find the
;;  best paramters.  (This is because the user probably didn't
;;  pick the exact minimum of Chi2.)
;;
;;  In addition, if extra errors are added to the data, this changes
;;  the location of the min Chi2.  (Presumably not much.)
;;  
;;  Do this N times, choosing initial points randomly around the
;;  Opti-Center, primarily because mpfit seems to miss the best fit.
;;
;; ------------------------------------------------------------

  N = 20

; 7 = 3 parms, 3 errors, Chi2
  parmgrid = fltarr(N, 7)

  For i = 0, N - 1 DO BEGIN

     ;; Seed Randon Numbers with SYSTIME(/SECONDS) + i (note -i wouldn't work)
     Rand = RandomU( Systime(/SECONDS) + i, 3)

     ThisSep = optiSep + (Rand[0] - .5)*SepRange 
     ThisAz  = optiAz  + (Rand[1] - .5)*AzRange
     ThisCrat = optiCrat + (Rand[2] - .5)*optiCrat*CratPRange 

     ;; ---- initial guess -----
     p0 = [ThisSep, ThisAz, ThisCrat ]

     ;; ---- constraints on the parameters ----
     pi = replicate({fixed:0, limited:[0,0], limits:[0,0], RelStep:.01}, 3)

       ;; Limit PA to within Range of Opti's

     pi[1].limited[0] = 1
     pi[1].limited[1] = 1
     pi[1].limits[0]  = OptiAz - AzRange
     pi[1].limits[1]  = OptiAz + AzRange

     pi[0].limited[0] = 1
     pi[0].limited[1] = 1
     pi[0].limits[0]  = OptiSep - SepRange
     pi[0].limits[1]  = OptiSep + SepRange 

     pi[2].limited[0] = 1
     pi[2].limits[0]  = 1.0

     ;; ---- Additional Arguments for 'binary_chi2' ----
     fa = {t3data:t3data}

     ; MPFit will return % Error in each parameter (if it's a good fit) 
     ; MPFit will return the Sum of Squares of values at best

     perror = fltarr(5)
     bestnorm = 0

     pout = mpfit('mp_binary_func_fast', p0, functargs=fa, parinfo=pi, perror=perror, bestnorm=bestnorm, status=status, errmsg=errmsg, quiet = 1)

     If status LE 0 then message, errmsg

     ;DOF = n_elements(t3data.t3phi) - 3
     DOF = 28 - 3
     NullDOF = 28
     errors = perror * SQRT(bestnorm/DOF)
     chi2 = bestnorm

     parmgrid[i,*] = [pout[0], errors[0], pout[1], errors[1], pout[2], errors[2], chi2]

     Print, i, ": ", Transpose( parmgrid[i,*] )
  EndFor

;; ------------------------------------------
;;   Final Chi2 and Parameter Errors
;;
;;   We get errors by noting that d^2(Chi2)/dparm^2 = 2/(err_parm)^2
;; ------------------------------------------

  FinalChi2     = Min( Parmgrid[*,6], Index )
  FinalSep     = Parmgrid[ Index, 0 ]
  FinalAz      = Parmgrid[ Index, 2 ]
  FinalCRat    = Parmgrid[ Index, 4 ]

  p0 = [ FinalSep, FinalAz, FinalCRat ]

  ; Errors 
  
  Errs = fltarr(3)
  percmove = .1
  For i=0, 3 - 1 DO BEGIN
     pplus = p0
     pplus[i] = pplus[i] * (1 + percmove)
     pneg  = p0
     pneg[i] = pneg[i] * (1-percmove)

     t3pMod = binary_t3data_fast( pplus, t3data=t3data )
     t3nMod = binary_t3data_fast( pneg , t3data=t3data )

     Chi2p = binary_chi2_fast( t3model=t3pmod, t3data=t3data )
     Chi2n = binary_chi2_fast( t3model=t3nmod, t3data=t3data )
     
     d2Chi2 = ( Chi2p - 2*FinalChi2 + Chi2n ) / ( percmove*p0[i] )^2
     d1Chi2 = ( Chi2p - Chi2n ) / (2*percmove*p0[i])
    
     Errs[i]   = Sqrt( (d1chi2/d2chi2)^2 + 2/d2chi2 )

  ENDFOR

  FinalSepErr = Errs[0]
  FinalAzErr  = Errs[1]
  FinalCRatErr= Errs[2]

  ;; Null Chi2
  ;;  Really should divide by 28, not 28-3 here.
  FinalNullChi2 = binary_chi2_fast( t3model=t3null, t3data=t3data )

  ;; Estimate of Confidence Level using F-Statistics
  ;;
  ;; The sum of two Chi2 variables (NOT reduced) is also a Chi2 variable with
  ;; degrees of freedom equal to the sum of the two.
  ;;
  ;;  DChi2 = NullChi2 - Chi2 is a Chi2 variable of degree 3
  ;;  F = DChi2 / Chi2_red
 

  DChi2 = ( FinalNullChi2 - FinalChi2 ) / 3.
  F = DChi2 / (FinalNullChi2 / NullDOF )

;; ----------------------------------------
;;
;;  This is a big hack.  Like, HHAACCKK
;;
;;  
;; -----------------------------------------

  return_info = {Sep: FinalSep, SepErr: FinalSepErr, Az: FinalAz, AzErr: FinalAzErr, Crat: FinalCrat, CratErr: FinalCratErr, Chi2: FinalChi2, NullChi2: FinalNullChi2, DataDate: "", AnalyseDate: "", Filter: "", Target: "", FileName: "" }

Median_Orig_Error = Median( t3data.t3phierr )

  Print, 'Median Original Error: ', Median_Orig_Error
  Print, ""
  Print, "Sep:  ", FinalSep, " +/- ", FinalSepErr
  Print, "Az:   ", FinalAz,  " +/- ", FinalAzErr
  Print, "CRat: ", FinalCrat," +/- ", FinalCRatErr
  Print, "Final Chi2: ", FinalChi2
  print, 'Null Chi2:  ', finalnullchi2
  Print, 'F-Statistic:  ', F, ' = ' , DChi2, ' / ', FinalNullChi2 / 28.
 

  
 t3best = binary_t3data_fast( [ optiSep, OptiAz, OptiCRat ], t3data=t3data )
 
 ;; Plot the best fit and residuals
 !p.multi = [ 0, 1, 2 ]
 s = sort( t3data.u3 )
 good = where( (t3data.t3phierr)[s] LT 100. )
 ploterr, good, ((t3data.t3phi)[s])[good], ((t3data.t3phierr)[s])[good], color=250, psym=5
 oplot, (t3best.t3phi)[s], color=65000

; plot, (t3data.t3phi / t3data.t3phierr )[good], color=250, /psym
; oplot, ((t3data.t3phi-t3best.t3phi)/t3data.t3phierr)[good], color=65000, /psym
 ploterr, ((t3data.t3phi-t3best.t3phi)[s])[good], ((t3data.t3phierr)[s])[good], color=250
 !p.multi = 0

;; -------------------------------
;;  Get Detection Confidence 
;;   & Companion Contrast Limits
;; -------------------------------

 FTestfile = file_dirname( datafile ) + "/" + file_basename( datafile, 'oifits' ) + "ftest.10000.idlvar"
 
 If( File_Test( FTestFile ) ) THEN BEGIN
    
    print, "Calculating Confidence of Detection..."
    restore, FTestFile

    FTest_parms = FTestInfo.Parmlimits
    Sim_CP = FTestInfo.Sim_CP
   
    emp_F_dist = FTestInfo.BestFitOutput[*,5]
    emp_F_dist = emp_F_dist[ Sort( emp_F_dist ) ]
    
    w = where( F GT emp_F_dist, psig )
    psig = 1.0 * psig / n_elements( emp_F_dist )

    print, "Significance of Detection: ", psig

    ;  Now, we want an achievable contrast limit for a given separation.
    ;   So, we average over the thetas.
    ;   Then, we extrapolate to what contrast limit will give us our desired
    ;   confidence limit.

    conf_limit = .995
    psig_cube = FTestInfo.psig_cube
    parm_cube = FTestInfo.parm_cube

    dummy = parm_cube[ 0, *, 0, 0 ] 
    FTestSep = reform( dummy, n_elements( dummy ) ) 
    n_sp = n_elements( FTestSep )
    dummy = parm_cube[ 2, 0, 0, * ]
    FTestCrat = reform( dummy, n_elements( dummy ) )
    n_cr = n_elements( FTEstCrat )
    psig_avg = total( psig_cube, 2 ) / n_cr
    Crat_Limits = fltarr( nSep )
    Mag_Limits  = fltarr( nSep )
    For i = 0, n_sp - 1 DO BEGIN
       ;; For Crat ~ 1.0, the signal is very small and so a Crat =
       ;; 1.01 is actually very difficult to measure.  So, psig_avg
       ;; peaks somewhere at Crat ~ a few.  We're interested in
       ;; only high crat detection probabilities, so cut out the low
       ;; ones.
       dummy = max( psig_avg[i,*], ind )
       
       Crat_lev = interpol( FTestCrat[ind:(n_cr-1)], psig_avg[i,ind:(n_cr-1)], conf_limit )

       ;ind_Crat is now the interpolated index of the Crat vector
       Crat_Limits[i] = Crat_lev
       Mag_limits[i]  = 2.5 * alog( Crat_lev ) / alog( 10. )
       ENDFOR
     Print, Transpose( [ [ FTestSep ], [ Crat_Limits ] ] )

 ENDIF


     

end
