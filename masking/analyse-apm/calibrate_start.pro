
pro main_event, event
  widget_control, event.top, get_uvalue=uvalue

  datafileinfo = uvalue.datafileinfo

  name_cube = datafileinfo.Directory + "/" + DataFileInfo.FileName
  root_dir  = uvalue.root_dir
  
  display_calibs_too = uvalue.display_calibs_too
  use_dwb_method     = uvalue.use_dwb_method
  skip_base_rejection = uvalue.skip_base_rejection

  calibrate_v2_cp, name_cube, root_dir, /reset, /average_src, display_calibs_too=display_calibs_too, use_dwb_method=use_dwb_method, skip_baseline_rejection=skip_base_rejection
END

pro calibs_too_event, event
  widget_control, event.top, get_uvalue=uvalue
  uvalue.display_calibs_too = event.select
  widget_control, event.top, set_uvalue=uvalue
END

pro use_dwb_method_event, event
  widget_control, event.top, get_uvalue=uvalue
  uvalue.use_dwb_method = event.select
  widget_control, event.top, set_uvalue=uvalue
END

pro skip_rejection_event, event
  widget_control, event.top, get_uvalue=uvalue
  uvalue.skip_base_rejection = event.select
  widget_control, event.top, set_uvalue=uvalue
END

pro CubeDropList_event, event
  ;; Take the value from the droplist and put it in the uvalue table
  widget_control, event.top, get_uvalue=uvalue
  uvalue.DataFileInfo = (uvalue.CubeInfo)[event.index]
  widget_control, event.top, set_uvalue=uvalue
END



pro Calibrate_Start


; Edittables
  root_dir  = '~/Research/idl-code/masking/'
  WORKINGDIR = "~/Research/redux/"
  ;WORKINGDIR = "~/Research/idl-code/SPF-SIM-Package/data_res/"
  CubeInfoFormatString = "cubeinfo*.idlvar"

  ;;
  ;; Grab the cubed data files, ready for calibration
  ;;
  
  CubeFiles = File_Search( WORKINGDIR, CubeInfoFormatString )
  BlankCubeStruct = { Target: '', Date: '', Filter: '', Directory: '', Filename: '' }

  CubeInfo = replicate( BlankCubeStruct, n_elements(CubeFiles) )

  For i = 0, n_Elements(CubeFiles) - 1 DO BEGIN
     FilePath = CubeFiles[i] 

     spl = strsplit( FilePath, "/", /Extract )

   ; Format: blah/Date/Target/Filter/cubeinfo*.idlvar
     n = n_elements( spl )
     Filter = spl[ n - 2 ]
     Target = spl[ n - 3 ]
     Date   = spl[ n - 4 ] 

     Struct = BlankCubeStruct
     Struct.Filter = Filter
     Struct.Target = Target
     Struct.Date   = Date
     Struct.Directory = FILE_DIRNAME( FilePath )
     Struct.FileName  = FILE_BASENAME( FilePath )

     CubeInfo[i] = Struct
  ENDFOR

  CubeInfo = Sort_Array_Of_Structs_By_Tags( CubeInfo, [ "Date", "Target", "Filter" ], 1)

  ;; ----------------------
  ;; Start Creating Widgets
  ;; ---------------------

  DataFileInfo = CubeInfo[0]
  UserValue = { CubeInfo: CubeInfo, DataFileInfo: DataFileInfo, root_dir: root_dir, display_calibs_too: 1, use_dwb_method: 1, Skip_Base_Rejection: 1 }

  ;; Base Widget
  baseID = widget_base( column = 1, title = 'Control Panel: Calibration', UValue=UserValue )

  ;; Droplist
  CubeDropListValues = replicate( '', n_elements( CubeInfo ) )
  For i = 0, n_elements(CubeInfo) - 1 DO CubeDropListValues[i] = (CubeInfo[i]).Date + " " + (CubeInfo[i]).Target + " " + (CubeInfo[i]).Filter
  
  CubedropListID = Widget_List(baseID, Value=CubedropListValues, event_pro='CubeDropList_event', ysize=20 )
  
  ;; Create the 'Display Calibs Too' Selection
  CalibsTooID = Widget_Base( baseID, row=1, /nonexclusive )
  CalibsTooOpt = Widget_Button( CalibsTooID, value="Display Calibs Too", event_pro='calibs_too_event' )
  Widget_Control, CalibsTooOpt, /Set_BUTTON

 ;; Create the 'USE DWB Method' Selection
  UseDWBID = Widget_Base( baseID, row=1, /nonexclusive )
  UseDWBOpt = Widget_Button( UseDWBID, value="Use DWB Method", event_pro='use_dwb_method_event' )
  Widget_Control, UseDWBOpt, /Set_BUTTON

 ;; Create the 'Skip Baseline Rejection' Selection
  SkipBaseID = Widget_Base( baseID, row=1, /nonexclusive )
  SkipBaseOpt = Widget_Button( SkipBaseID, value="Skip Baseline Rejection", event_pro='skip_rejection_event' )
  Widget_Control, SkipBaseOpt, /Set_BUTTON

  ;; Go Bots Go!  Refresh!
  RefreshID = widget_button( baseID, value="Calibrate!" )

  Widget_Control, baseID, /realize
  xmanager, 'Calibrate_Start', baseID, event_handler = 'main_event'

END

