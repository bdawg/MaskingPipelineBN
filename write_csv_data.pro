;*****************************************************************************************************
;+
; NAME:
;       WRITE_CSV_DATA
;
; PURPOSE:
;
;       This procedure allows the user to write a 2D array in a format
;       suitable for inclusion in a spreadsheet. By default it produces
;       a comma separated ASCII data file, but any other delimiter can
;       also be used.
;
; AUTHOR:
;
;        FANNING SOFTWARE CONSULTING
;        1645 Sheely Drive
;        Fort Collins
;        CO 80526 USA
;        Phone: 970-221-0438
;        E-mail: davidf@dfanning.com
;
; CATEGORY:
;
;       File I/O.
;
; CALLING SEQUENCE:
;
;       Write_CSV_Data, theData
;
; ARGUMENTS:
;
;        data:            A 2D array of data to write in spreadsheet format.
;
; KEYWORDS:
;
;        COLUMNHEADERS:   A vector of column headers. Must be the same size as the X dimension
;                         of the data array. This argument is optional.
;
;        DELIMITER:       The type of delimiter to use between array values. By default, a comma (",").
;                         If you wanted to use a tab, for example, you would set DELIMITER=String(9B).
;                         The value must be a string scalar.
;
;        FILENAME:        The name of the output file. Will be requested of the user if not passed in.
;
;        WIDTH:           The width of the output line. The default is 1600 characters. The width
;                         must be wide enough to accommodate a row of the data array, written as
;                         a string vector (with deliminators).
;
; MODIFICATION_HISTORY:
;
;       Written by: David W. Fanning, 12th March 2003.
;       DELIMETER keyword added by Ben Tupper, 7 January 2004.
;-
;*****************************************************************************************************
PRO Write_CSV_Data, $
   data, $
   ColumnHeaders=columnHeaders, $
   Filename=filename, $
   Width=width, $
   Delimiter = delimiter

   ; Return to caller on error.

On_Error, 1

   ; Must have data to write.

IF N_Elements(data) EQ 0 THEN BEGIN
   Message, 'DATA parameter is required. Returning...', /Informational
   RETURN
ENDIF

   ; Get the number of dimensions of the data. Data must be 2D.

ndims = Size(data, /N_Dimensions)
IF ndims NE 2 THEN BEGIN
   Message, 'DATA must be 2D. Returning...', /Informational
   RETURN
ENDIF

   ; Get the sizes of the 2D data.

s = Size(data, /Dimensions)
xsize = s[0]
ysize = s[1]

   ; Check for column header vector.

IF N_Elements(columnHeaders) NE 0 THEN BEGIN
    length = N_Elements(columnHeaders)
    IF length NE xsize THEN BEGIN
      Message, 'Column Header vector wrong size. Returning...', /Informational
      RETURN
    ENDIF
ENDIF

   ; Get a filename if you need one.

IF N_Elements(filename) EQ 0 THEN filename = Dialog_Pickfile(/Write, File='data.csv')
IF filename EQ "" THEN RETURN

   ; Check keywords.

IF N_Elements(width) EQ 0 THEN width = 1600

   ; Need comma separated data file.
   ; or some variant of the comma
IF n_elements(delimiter) NE 0 THEN comma = delimiter[0] ELSE comma = ","

   ; Open the data file for writing.

OpenW, lun, filename, /Get_Lun, Width=width

   ; Write the column header row, if available.

IF N_Elements(columnHeaders) NE 0 THEN BEGIN

      ; Make sure these are strings.

   sColumns = StrTrim(columnHeaders, 2)

      ; Add a comma to each value except the last one.

   sColumns[0:xsize-2] = sColumns[0:xsize-2] + comma

      ; Write the headers to the file.

   PrintF, lun, sColumns

ENDIF

   ; Write the data to the file.

sData = StrTrim(data,2)
sData[0:xsize-2, *] = sData[0:xsize-2, *] + comma
PrintF, lun, sData

   ; Close the file.

Free_Lun, lun
END