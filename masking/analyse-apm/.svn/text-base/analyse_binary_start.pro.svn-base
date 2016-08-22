;; -----------------------------------------------------------------------------
;;                           FIT A BINARY SYSTEM
;; -----------------------------------------------------------------------------

;; -------------------- read the data oifits --------------------

;datafile = '~/Research/redux/070829/G_244/H/G_244-049_F_1580.oifits'
;datafile = '~/Research/redux/060213/G_244/H/G_244-049_F_0254.oifits'
;datafile = '~/Research/redux/061010/G_244/H/G_244-049_F_0042.oifits'
;datafile = '~/Research/redux/061010/G_244/K_short/G_244-049_F_0062.oifits'
;datafile = '~/Research/redux/061010/G_244/K_short/G_244-049_F_0142.oifits'
;datafile = '~/Research/redux/070731/G_212/H/G_212-057_F_1762.oifits'
;datafile = '~/Research/redux/070926/G_212/H/G_212-057_F_0508.oifits'
;datafile = '~/Research/redux/070829/G_212/H/2M2146+4638_F_1036.oifits'
;datafile = '~/Research/redux/070829/GJ_101/H/GJ_101_F_1780.oifits'
;datafile = '~/Research/redux/060214/GJ_101/H/GJ_101_F_0223.oifits'

pro Analyse_event, event
  widget_control, event.top, get_uvalue=uvalue
  
  datafileinfo = uvalue.datafileinfo
  datafile = datafileinfo.directory + "/" + datafileinfo.filename

  parmlimits = uvalue.parmlimits

  analyse_binary, datafile, parmlimits, return_info

  ;; This is the information returned from analyse_binary in a handy
  ;; package.

  return_info.DataDate = datafileinfo.Date
  return_info.Filter = datafileinfo.Filter
  return_info.Target   = datafileinfo.Target
  return_info.FileName = datafileinfo.Filename
  return_info.AnalyseDate = systime( 0 )
END

pro ParmTable_event, event
  ;; Just get the table values and store them in the top uvalue
  widget_control, event.id, get_value=table
  widget_control, event.top, get_uvalue=uvalue
  uvalue.parmlimits = table
  widget_control, event.top, set_uvalue=uvalue
END

pro OIFitsDropList_event, event
  ;; Take the value from the droplist and put it in the uvalue table
  widget_control, event.top, get_uvalue=uvalue
  uvalue.DataFileInfo = (uvalue.OIFitsInfo)[event.index]
  widget_control, event.top, set_uvalue=uvalue
END

pro analyse_binary_start
  datafile = "Blah"
  WORKINGDIR = "~/Research/redux/"
  ;;WORKINGDIR = "~/Research/idl-code/SPF-SIM-Package/data_res/"
  ;;
  ;; Before creating widgets, grab the calibrated data, ready for analysis.
  ;;
  
  OIFitsFiles = File_Search( WORKINGDIR, "*.oifits" )
  BlankCubeStruct = { Target: '', Date: '', Filter: '', Directory: '', Filename: '', HasDWB: '', HasFTest: '' }

  FirstIndexFlag = 0
  FOR i = 0, n_elements( OIFitsFiles) - 1 DO BEGIN
     FilePath = OIFitsFiles[i] 

     ;;If( Strpos( filepath, "calib_0809" ) NE -1 ) THEN Continue

     spl = strsplit( FilePath, "/", /Extract )
     
     ; Format: blah/Date/Target/Filter/[Target_*.oifits]
     n = n_elements( spl )
     Filter = spl[ n - 2 ]
     Target = spl[ n - 3 ]
     Date   = spl[ n - 4 ]
     FileName = FILE_BASENAME( FilePath) 
     
     ;; There are other oifits files around.  Visibility and CP Data
     ;; file names are formatted such that the OIFITS
     ;; file name starts with the Target Name.  Then, we're gold.

     ;FormatTarget = Format_Target_Name( Target )
     ;If( StrPos( FileName, FormatTarget ) EQ -1 ) THEN Continue

     If( StrPos( Filename, ".vis." ) NE -1 ) Then Continue

     Struct = BlankCubeStruct
     Struct.Filter = Filter
     Struct.Target = Target
     Struct.Date   = Date
     Struct.Directory = FILE_DIRNAME( FilePath )
     Struct.FileName  = FileName

     If( StrPos( Filename, ".dwb." ) NE -1 ) THEN BEGIN
        Struct.HasDWB = 'DWB'
        If( n_elements( FILE_SEARCH( Struct.Directory, "*.dwb.ftest.10000.*" ) ) NE 0 ) THEN Struct.HASFTest = '*'
     ENDIF

     If( StrPos( Filename, ".xxx.nmf." ) NE -1 ) THEN BEGIN
        Struct.HasDWB = 'NMF'
        If( n_elements( FILE_SEARCH( Struct.Directory, "*.xxx.nmf.ftest.10000.*" ) ) NE 0 ) THEN Struct.HASFTest = '*'
     ENDIF

     ;; Skip Non-DWB
     ;;If( Struct.HasDWB EQ '' ) THEN Continue
     ;;If( Struct.HasDWB NE '' ) THEN Continue

     If( FirstIndexFlag EQ 0 ) THEN BEGIN
        OIFitsInfo = [ Struct ]
        FirstIndexFlag = 1
     ENDIF ELSE OIFitsInfo = [ Struct, OIFitsInfo ] 
  ENDFOR

  OIFitsInfo = Sort_Array_Of_Structs_By_Tags( OIFitsInfo, [ "Date", "Target", "Filter" ], 1)

  ;;
  ;; This is the structure of the opaque parmlimits 
  ;; Index: 0 = sep, 1 = Az, 2 = Crat
  ;;

  parmlimits = [ { min: 45., max: 250., n: 15 }, { min: 0., max: 380., n: 20 }, { min: 0., max: 3.5, n: 15} ] 
  parmlimits = [ { min: 45., max: 450., n: 20 }, { min: 0., max: 380., n: 20 }, { min: 0., max: 1.9, n: 20} ]
  
  ;; Assume the initil datafile is the first one on the list

  DataFileInfo = OIFitsInfo[0]

  ;; 
  ;; Start Creating Widgets
  ;; 

  UserValue = { OIFitsInfo: OIFitsInfo, DataFileInfo: DataFileInfo, ParmLimits: ParmLimits  }

  ;; The Base Widget
  BaseID = widget_base( column=2, UValue = UserValue, title = "Control Panel: Analyse" )

  ;; The ParmLimit Table Widget
  TableID = widget_table( BaseID, value = parmlimits, uvalue=UserValue, /editable , column_labels = Tag_Names(parmlimits), row_labels = [ "Sep: ", "Az:  ", "CRat:" ], event_pro='parmtable_event' )

  ;; The OIFits DropList Widget
  OIFitsDropListValues = replicate( '', n_elements( OIFitsInfo ) )
  For i = 0, n_elements(OIFitsInfo) - 1 DO OIFitsDropListValues[i] = (OIFitsInfo[i]).Date + " " + (OIFitsInfo[i]).Target + " " + (OIFitsInfo[i]).Filter + " " + (OIFitsInfo[i]).HasDWB + " " + (OiFitsInfo[i]).HasFTEst

  OIFitsdropListID = Widget_List( baseID, Value=OIFitsDropListValues, event_pro='OIFitsDropList_Event', ysize=20 , /multiple )
  
  ;; Go Bots Go!  Analyse!
  RefreshID = widget_button( baseID, value="Analyse" )

  widget_control, baseID, /realize

  ;; Use a Window As Large As Possible
  device, get_screen_size=screen_size
  window, /free, xsize=screen_size[0]*.65, ysize=screen_size[1]*.9 
 
  xmanager, 'Analyse_Binary_Start', baseID, event_handler = 'Analyse_event'

END
