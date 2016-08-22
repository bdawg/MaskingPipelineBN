pro cron_emp_f_dist, WORKINGDIR = WORKINGDIR

  If( NOT Keyword_Set(WORKINGDIR) ) THEN WORKINGDIR = "~/Research/redux/"

  WorkingString = "dwb.is.working"
  WorkingFileText = "This file is a placeholder to notify parallel processes that work is already being done here.  If for some reason you need to track down its owner, a disasterous error has probably occurred, so I wont make it any easier to find me than it already is."

  ;; Standard Parameter Limits (for "Non-Near" Binaries)
  MDwarf_parmlimits = [ { min: 43., max: 240., n: 15 }, { min: 0., max: 380., n: 15 }, { min: 0., max: 2.8, n: 15} ]
  LDwarf_parmlimits = [ { min: 43., max: 240., n: 15 }, { min: 0., max: 380., n: 15 }, { min: 0., max: 1.4, n: 15} ]
  

  NSim = 10000
  
;; -------------------------- Start of Code ----------------------
 
  ;; First, find all targets that have dwb.extra.idlvar
  ;;  and also don't have dwb.ftest.10000.idlvar
  
  NSimString = String( NSim, format='(I5.5)' )
  ExtraFiles  = File_Search( WORKINGDIR, "*.dwb.extra.idlvar" )
  ExtraDIR    = FILE_DIRNAME( ExtraFiles )
  OiFitsFiles = File_Search( WORKINGDIR, "*.dwb.oifits" )
  OiFitsDIR   = FILE_DIRNAME( OiFitsFiles )

  FTestFiles  = File_Search( WORKINGDIR, "*.dwb.ftest." + NSimString + ".*" )
  FTestDIR    = FILE_DIRNAME( FTestFiles )

  ;; Hack Out Anything But L Dwarfs For Monday.
  S = WHERE( -1 NE Strpos( ExtraDIR, '2M' ) )
  ExtraFiles = ExtraFiles[S]
  ExtraDIR   = ExtraDir[S]

  ;; Using Craig Markwardt CMSET_OP.pro code (see his website). 
  ;; And the directories that CanDo but not already done
  CanDO   = CMSET_OP( ExtraDIR, 'AND', OiFitsDIR )
  ToDoDIR = CMSET_OP( CanDo, 'AND', /NOT2, FTestDIR )

  ;;  Now, go through each one and do it, as long as another
  ;;  process isn't already working on it.  When a process
  ;;  starts working on it, it creates a small .working file.

  ErrorStack = [ "The following Errors were found: " ]
  TestsRun = 0
  FOR n_test = 0, n_elements( ToDoDIR ) - 1 DO BEGIN

     If( TestsRun EQ 10 ) THEN BEGIN
        Msg = "Ran 10 F-Tests This Run.  Calling is Quits.  Congrats.  Grab a drink."
        ErrorStack = [ ErrorStack, MSg ]
        break
     ENDIF

     FilePath = ToDoDIR[n_test] 

     ; Check to see if a process is already working here
     IsWorking = File_Search( FilePath, WorkingString, count=cw )
     If( cw NE 0 ) THEN BEGIN
        Msg = "WARNING: The working file '" + WORKINGString + "' is found in " + FilePath + ".  Skipping."
        print, Msg
        ErrorStack = [ ErrorStack, Msg ]
        Continue
     ENDIF

     ; Eh, sue me.
     ExtraFile  = File_Search( FilePath, "*.dwb.extra.idlvar", count=c1 )
     OiFitsFile = File_Search( FilePath, "*.dwb.oifits", count=c2 )
     
     If( (c2 GT 1) OR (c1 GT 1) ) THEN BEGIN
        Msg = "ERROR: More than one ExtraFile or dwb.oifits file found in " + FilePath + ".  Skipping."
        print, Msg
        ErrorStack = [ ErrorStack, Msg ]
        Continue
     ENDIF

     ;; ------------------------------
     ;;  Start the long haul
     ;; -----------------------------

     ;; Create Placeholder File
     Command = "echo '" + WorkingFileText + "' > " + FilePath + "/" + WorkingString
     Print, "Running UNIX Command: ", Command
     Spawn, Command

     print, replicate( "****", 1, 3 ), "Loading DWB.OIFITS File for Directory: ", FilePath, replicate( "****", 1, 3 )

     ; ExtraFile contains variable Bispect and "bad_bispect" and TargetInfo
     ;; HACK here - Actually, a few lines for backwards compatibility 
     TargetInfo = { TargetType: 'm dwarf' } ;; Older versions of calibrated data do not have this tag, default to m dwarf.
     cp_cov = 0.0
     
     restore, ExtraFile
     cp_cov = full_stats.science_stats.cp_cov

     If( TargetInfo.TargetType EQ 'l dwarf' ) THEN parmlimits = LDwarf_parmlimits ELSE parmlimits = MDwarf_parmlimits
     Print, "Star Type: ", TargetInfo.TargetType

     extract_t3data, file=oifitsfile, t3data
      
     ;; ------------------------------------------
     ;;
     ;;  Run emp_f_dist.pro on the DWB.OIFITS Data
     ;;
     ;; -----------------------------------------

     print, replicate( "****", 1, 3 ), "Calculating F-Test Distribution for Directory: ", FilePath, replicate( "****", 1, 3 )
     
     ; Go
     emp_F_dist, t3data, parmlimits, cp_cov, ret_struct, nSim=NSim

     ;; ---------------------------------------
     ;;
     ;; Calculate the achievable Contrast Limits
     ;;  from the ftest data, save it in the ftest file
     ;;
     ;; ---------------------------------------

     emp_f_dist = ret_Struct.BestFitOutput[*,5]
     emp_f_dist = emp_f_dist[ Sort( emp_f_dist ) ]
     Sim_CP = ret_struct.Sim_CP
     NullChi2 = (ret_struct.BestFitOutput)[*,3]

     print, replicate( "****", 1, 3 ), "Calculating Achievable Contrast Limits for Files in Directory: ", FilePath, replicate( "****", 1, 3 )
     
     emp_crat_limits, t3data, parmlimits, Sim_CP, emp_F_dist, nullchi2, ret_struct_crat

     FTestInfo = { parmlimits: parmlimits, cp_cov: cp_cov, t3data: t3data, nSim: nSim, BestFitOutput: ret_Struct.BestFitOutput, Sim_CP: ret_Struct.Sim_CP, psig_cube: ret_struct_crat.psig_cube, parm_cube: ret_struct_crat.parm_cube, version: "Version 1.0" }

     FTestfile = FilePath + "/" + file_basename( oifitsfile, '.oifits' ) + ".ftest." + NSimString + ".idlvar"
     
     ; Save Data
     save, FTestInfo, filename=FTestFile

     ;; --------------------------
     ;;  Clean Up
     ;; --------------------------

     Command = "rm " + FilePath + "/" + WorkingString
     Print, "Running Command: ", Command
     Spawn, Command

     Msg = "Successfully calculated F-Distribution (N=" + NSimString + ") for " + FilePath
     ErrorStack = [ ErrorStack, Msg ]
     TestsRun++
  ENDFOR

  Print, "Let's call it a night."
  Print, ErrorStack

END
