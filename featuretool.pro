;+
; NAME:
;     FEATURETOOL
;
; PURPOSE:
;     Graphical interface for setting parameters for FEATURE.
;
; CATEGORY:
;     Image processing, particle tracking
;
; CALLING SEQUENCE:
;     IDL> featuretool
;
; OPTIONAL INPUTS:
;     a: grey-scale image
;  
; SIDE EFFECTS:
;     Opens a window.
;
; RESTRICTIONS:
;     Only works properly on grayscale images.
;
; PROCEDURE:
;     Calls BPASS, FEATURE, and DEINTERLACE according to user inputs.
;
; MODIFICATION HISTORY:
; 12/15/2005: Written by David G. Grier, New York University
; 05/19/2006: DGG. Added deinterlace and invert.  Reorganized.
;    Defaults to pgm/ppm.  Understands more file formats.
; 05/23/2006: DGG. Added region of interest.  Overlay image with
;    region of interest and extent of BPASS.  Save features and
;    save procedure.  More comments.
; 06/05/2006: DGG. Corrected region-of-interest code. Reorganized
;    automatically generated procedure code.
; 06/13/2006: DGG. Improved automatically generated procedure code.
;    Added "Save Image" feature.  Added IDLSNAP feature.
;    Made non-blocking.  
;    SEPARATION must exceed DIAMETER by at least 1.
;    Fixed limits when image size changes.
; 03/29/2010: DGG. Fixed DIALOG_PICKFILE bug.  Added rudimentary
;    insructions.  Corrected quotation marks throughout.  Formatting.
; 06/09/2010: DGG. Support for input image on command line.  IDLSNAP
;    desensitized on systems lacking idlsnap support.  Fixed file
;    specification in emitted code.  Fixed crash when no features are found.
; 06/10/2010: DGG. Added COMPILE_OPT statements.  Corrected IDLSNAP
;    behavior when no video signals are detected.  Separate and format
;    help text and metaprogramming code.
;
; Copyright (c) 2005-2010 David G. Grier
;-

;;;
;;; FEATURETOOL_INSTRUCTIONS
;;;
;;; Open a dialog box with instructions
;;;
function featuretool_instructions

COMPILE_OPT IDL2, HIDDEN

instructions = [$
'1. Use the dialog to select and open an image file.', $
'   Otherwise, File->IDLSNAP to acquire an image from a camera, if implemented on your system.', $
'2. Adjust XMIN, XMAX, YMIN, and YMAX to define region of interest, if desired.', $
'3. Select Deinterlace and Invert, if needed.', $
'4. Select BPASS tab.', $
'4.a. Adjust EXTENT until filtered image looks clean, not blurred.', $
'4.b. Adjust NOISE to clean up noisy pixels.  This usually is not necessary', $
'5. Select IMAGE tab.', $
'5.a. Adjust DIAMETER to size of objects, in pixels.', $
'5.b. Adjust SEPARATION to minimum separation between objects.', $
'5.c. Adjust MIN to eliminate false local maxima.', $
'At this point the identified features should look good.', $
'6. Select MASSCUT tab and make sure that all identified features form a compact cloud of data points.', $
'   If not, adjust MASSCUT to eliminate bright outliers.', $
'7. Select DIAMETER tab to make sure that the histogram of particle offsets is reasonably flat.', $
'   A dip in the middle indicates pixel bias; increase DIAMETER until the distribution is flat.', $
'7 File->Save Feature: Save information about the features identified in this image.', $
'8 File->Save Routine: Save an IDL procedure that implements the present analysis.']
res = DIALOG_MESSAGE(instructions, /INF, TITLE='Feature Tool: Instructions')

return, res
end

;;;
;;; FEATURETOOL_GETVALUE
;;;
;;; Get value from a widget, making sure that the value
;;; remains within bounds
;;;
function featuretool_getvalue, wid, min, max

COMPILE_OPT IDL2, HIDDEN

  WIDGET_CONTROL, wid, GET_VALUE = val
  val = val < max > min
  WIDGET_CONTROL, wid, SET_VALUE = val
  return, val
end

;;;
;;; FEATURETOOL_UPDATE_ROI
;;;
;;; Update widget info regarding the region of interest
;;;
pro featuretool_update_roi, s, a

COMPILE_OPT IDL2, HIDDEN

  sz = size(a, /dimensions)
  ; reset ROI limits if image size has changed
  if s.p.w ne sz[0] or s.p.h ne sz[1] then begin
     s.p.w = sz[0]
     s.p.h = sz[1]
     s.p.xmin = 0
     s.p.ymin = 0
     s.p.xmax = s.p.w - 1
     s.p.ymax = s.p.h - 1
  endif else begin
     s.p.xmax = s.p.xmax < (s.p.w - 1)
     s.p.ymax = s.p.ymax < (s.p.h - 1)
     s.p.xmin = s.p.xmin < (s.p.xmax - 1)
     s.p.ymin = s.p.ymin < (s.p.ymax - 1)
  endelse
  WIDGET_CONTROL, s.w.wxmin, SET_VALUE = s.p.xmin
  WIDGET_CONTROL, s.w.wxmax, SET_VALUE = s.p.xmax
  WIDGET_CONTROL, s.w.wymin, SET_VALUE = s.p.ymin
  WIDGET_CONTROL, s.w.wymax, SET_VALUE = s.p.ymax
end

;;;
;;; FEATURETOOL_PROCESSIMAGE
;;;
;;; Apply all selected operations to the current image
;;; and update widgets with the results
;;;
function featuretool_processimage, s

COMPILE_OPT IDL2, HIDDEN

                               ; get a clean copy of the image
   a = read_image(s.file)      ; why?  This already should be clean
   
   ; preprocess image
   if s.p.deinterlace then a = deinterlace(a)
   if s.p.invert then a = 255b - a

   ; update region of interest
   featuretool_update_roi, s, a
    
   ; BPASS on region of interest. Do it this way around to avoid
   ; artifacts due to bad regions in original image
   b = bpass(a[s.p.xmin:s.p.xmax,s.p.ymin:s.p.ymax], $
             s.p.noise, s.p.extent)
                                ; at last!  Find the features.
   f = feature(b, s.p.diameter, s.p.separation, $
               masscut = s.p.masscut, min = s.p.min, $
               /iterate, /quiet, lmax=lmax, count=count)
   s.p.nlmax = lmax             ; number of local maxima
   s.p.nfeatures = count        ; number of qualified features

   if s.p.nfeatures gt 0 then begin
      ; index features to original image, not clipped image
      f[0,*] += s.p.xmin
      f[1,*] += s.p.ymin
   endif

   ; overlay features on image
   WIDGET_CONTROL, s.w.wimage, GET_VALUE = ndx
   wset, ndx
   plotimage, bytscl(a), /PRESERVE_ASPECT, xmargin = [5, 1], ymargin = [2, 1]
   if s.p.nfeatures gt 0 then $
      oplot, f[0, *], f[1, *], psym = circ()
   ; overlay clipping region of interest on image
   roi = [[s.p.xmin,s.p.ymin], $
          [s.p.xmin,s.p.ymax], $
          [s.p.xmax,s.p.ymax], $
          [s.p.xmax,s.p.ymin], $
          [s.p.xmin,s.p.ymin]]
   plots,roi
   ; overlay border corresponding to BPASS extent
   roi = [[s.p.xmin+s.p.extent,s.p.ymin+s.p.extent], $
          [s.p.xmin+s.p.extent,s.p.ymax-s.p.extent], $
          [s.p.xmax-s.p.extent,s.p.ymax-s.p.extent], $
          [s.p.xmax-s.p.extent,s.p.ymin+s.p.extent], $
          [s.p.xmin+s.p.extent,s.p.ymin+s.p.extent]]
   plots,roi,linestyle=2

   ; bpass-filtered image
   WIDGET_CONTROL, s.w.wbpass, GET_VALUE = ndx
   wset, ndx
   plotimage, bytscl(b), /PRESERVE_ASPECT, xmargin = [5, 1], ymargin = [2, 1]

   ; radius of gyration versus integrated brightness
   WIDGET_CONTROL, s.w.wmasscut, GET_VALUE = ndx
   wset, ndx
   if s.p.nfeatures gt 0 then $
      plot, f[2, *], f[3, *], psym = 3, $
            xtitle = 'm', ytitle = 'R_g', charsize = 1.5 $
   else $
      plot, [0, 1], /nodata, $
            xtitle = 'm', ytitle = 'R_g', charsize = 1.5
    
   ; distribution of fractional pixel positions -- should be flat
   WIDGET_CONTROL, s.w.wdiameter, GET_VALUE = ndx
   wset, ndx
   binsize = 0.05
   hist = s.p.hist
   nbins = n_elements(hist)
   if s.p.nfeatures gt 0 then $
      hist = histogram(f[0:1, *] mod 1, $
                       min = 0, max = 1, nbins = nbins, $
                       locations = x, input = hist)
   plot, x, hist, $
         xtitle = 'f_j', ytitle = 'N(f_j) df_j', charsize=1.5
   s.p.hist = hist

   ; number of local maxima
   WIDGET_CONTROL, s.w.wlmax, SET_VALUE = s.p.nlmax
    
   ; number of features
   WIDGET_CONTROL, s.w.wfeat, SET_VALUE = s.p.nfeatures

   return, f
end

;;;
;;; FEATURETOOL_SAVEROUTINE
;;;
;;; Metaprogramming: Save an IDL routine that analyzes images
;;; based on the present settings.
;;;
pro featuretool_saveroutine, s

COMPILE_OPT IDL2, HIDDEN

  fn = DIALOG_PICKFILE(DEFAULT_EXTENSION = 'pro', FILTER = ['*.pro'], $
                       /WRITE, /OVERWRITE_PROMPT, TITLE = 'Save Procedure')
  if strlen(fn) lt 1 then return

  nl = (!d.name eq 'WIN') ? string([13b, 10b]) : string(10b) ; newline

  openw, lun, fn, /get_lun
  printf, lun, $
 ';;; routine generated automatically by FEATURETOOL.PRO' + nl + $
 ';;; ' + systime()                                       + nl + nl + $
 ';;; Specify the image files to be processed:'           + nl + $
 ';;;    Example: filespec = "junk/*/*/*.pgm"'            + nl + $
 'pro '+file_basename(fn, '.pro')+', filespec'            + nl + nl + $
 'if n_params() eq 0 then filespec = "'+s.file+'"'        + nl

  if ((s.p.xmax - s.p.xmin + 1) ne s.p.w) or $
     ((s.p.ymax - s.p.ymin + 1) ne s.p.h) then begin
     printf, lun, $
 ';;; Region of interest'                                 + nl + $
 'xmin = ' + strtrim(s.p.xmin, 2)                         + nl + $
 'xmax = ' + strtrim(s.p.xmax, 2)                         + nl + $
 'ymin = ' + strtrim(s.p.ymin, 2)                         + nl + $
 'ymax = ' + strtrim(s.p.ymax, 2)                         + nl
  endif
  
  printf, lun, $
 ';;; Parameters for BPASS'                               + nl + $
 'noise = ' + strtrim(s.p.noise, 2)                       + nl + $
 'extent = ' + strtrim(s.p.extent, 2)                     + nl

  printf, lun, $
 ';;; Parameters for FEATURE'                             + nl + $
 'diameter = ' + strtrim(s.p.diameter, 2)                 + nl + $
 'separation = ' + strtrim(s.p.separation, 2)             + nl + $
 'min = ' + strtrim(s.p.min, 2)                           + nl + $
 'masscut = ' + strtrim(s.p.masscut, 2)                   + nl

  printf, lun, $
 ';;; Process file ...'                                   + nl + $
 'f = file_search(filespec, count = nfiles)'              + nl + $
 'for n = 0, nfiles-1 do begin'                           + nl + $
 '   a = read_image(f[n])'

  if ((s.p.xmax - s.p.xmin + 1) ne s.p.w) or $
     ((s.p.ymax - s.p.ymin + 1) ne s.p.h) then $
        printf, lun, $
 '   a = a[xmin:xmax,ymin:ymax]'

  if s.p.deinterlace then $
     printf, lun, $
 '   a = deinterlace(a)'

  if s.p.invert then $
     printf, lun, $
 '   a = 255b - a'

  printf, lun, $
 '   b = bpass(a, noise, extent)'                         + nl + $
 '   p = feature(b, diameter, separation, $'              + nl + $
 '               min=min, masscut=masscut, /iterate, /quiet)'  + nl + $
 '   write_gdf, p, f[n]+".gdf"'                           + nl + $
 'endfor'                                                 + nl + $
 'end'
  close, lun
  free_lun, lun
end

;;;
;;; FEATURETOOL_EVENT
;;;
;;; Process the event loop
;;;
pro featuretool_event, ev

COMPILE_OPT IDL2, HIDDEN

; The base widget's uvalue is the structure containing all the
; information about the program's state.
WIDGET_CONTROL, ev.TOP,  GET_UVALUE = s

; Clicking on a tab is handled automatically.  Nothing to do, so
; move along.
if (ev.ID eq s.w.wtab) then return

; Other widgets can give us something to do. 
WIDGET_CONTROL, ev.ID, GET_UVALUE = uval

doimage = 0                     ; Assume we're not reading a new image ...
doclearhist = 1                 ; ... so we should clear the feature histogram.
savefeatures = 0                ; Write features to a file.
CASE uval of
   'FILE': begin
      doclearhist = 0           ; Getting a new image
      CASE ev.DONE of
         0: begin               ; Still entering filename
            s.file = ev.VALUE
            WIDGET_CONTROL, ev.TOP, SET_UVALUE=s
         end
         1: begin               ; We have a complete file name ...
            s.havefile = 0      ; ... but we don't yet know it's good
            doimage = 0
            if query_image(s.file) then begin
               s.havefile = 1
               doimage = 1
            endif $
            else begin
               txt = [s.file, 'is not a recognized image file']
               res = DIALOG_MESSAGE(txt, /INF, $
                                    TITLE = 'Unrecognized Image Format')
               s.file = ''
            endelse
            WIDGET_CONTROL, ev.TOP, SET_UVALUE=s
         end
         2:
      endcase
   end

   'XMIN': begin
      s.p.xmin = featuretool_getvalue(ev.ID, $
                                      0, s.p.xmax - 1)
      doimage = 1
   end

   'XMAX': begin
      s.p.xmax = featuretool_getvalue(ev.ID, $
                                      s.p.xmin + 1, s.p.w - 1)
      doimage = 1
   end
   'YMIN': begin
      s.p.ymin = featuretool_getvalue(ev.ID, $
                                      0, s.p.ymax - 1)
      doimage = 1
   end
   'YMAX': begin
      s.p.ymax = featuretool_getvalue(ev.ID, $
                                      s.p.ymin + 1, s.p.h - 1)
      doimage = 1
   end
   
   'BPPREPROCESS': begin
      CASE ev.value of
         'BPDEINTERLACE': begin
            s.p.deinterlace = ev.select
            doimage = 1
         end
         'BPINVERT': begin
            s.p.invert = ev.select
            doimage = 1
         end
      endcase
   end
   'BPNOISE': begin
      s.p.noise = featuretool_getvalue(ev.ID, 0, 5)
      doimage = 1
   end
   'BPEXTENT': begin
      s.p.extent = featuretool_getvalue(ev.ID, 1, 25)
      doimage = 1
   end
   'FDIAMETER': begin
      s.p.diameter = featuretool_getvalue(ev.ID, 1, 25)
      doimage = 1
   end
   'FSEPARATION': begin
      s.p.separation = featuretool_getvalue(ev.ID, 2, 100)
      doimage = 1
   end
   'FMASSCUT': begin
      s.p.masscut = featuretool_getvalue(ev.ID, 0, 5000)
      doimage = 1
   end
   'FMIN': begin
      s.p.min = featuretool_getvalue(ev.ID, 1, 254)
      doimage = 1
   end
   
   'SAVEIMAGE': begin
      if s.havefile then begin
         fn = DIALOG_WRITE_IMAGE(read_image(s.file),/WARN_EXIST)
      endif else begin
         msg = ['No image to save', 'Open an image file']
         res = DIALOG_MESSAGE(msg, /INF)
      endelse 
   end

   'SAVEFEATURES': begin
      if s.havefile then begin
         fn = DIALOG_PICKFILE(/WRITE, /OVERWRITE_PROMPT, $
                              TITLE = 'Save Features', $
                              FILTER = ['*.gdf'], DEFAULT_EXTENSION = 'gdf')
         savefeatures = strlen(fn) GE 1
         doimage = savefeatures
      endif else begin
         msg = ['No features to save', 'Open an image file']
         res = DIALOG_MESSAGE(msg, /INF)
      endelse
   end

   'SAVEPRO': begin
      featuretool_saveroutine, s
   end
   
   'IDLSNAP': begin
      clearhist = 0
      a = idlsnap(/quiet)
      if n_elements(a) gt 1 then begin
         s.file = filepath('idlsnap.pgm',/TMP)
         write_ppm, s.file, a
         s.havefile = 1
         doimage = 1
      endif else begin
         s.havefile = 0
         doimage = 0
         instructions = 'Cannot snap photo.  Check video source'
         res = DIALOG_MESSAGE(instructions, /INF, TITLE = 'IDLSNAP ERROR')
      endelse
   end

   'HELPINSTRUCTIONS': begin
      res =  featuretool_instructions()
   end
   
   'HELPABOUT': begin
      about = ['FEATURETOOL v. 1.1', $
               'Written by David G. Grier', $
               'New York University']
      res = DIALOG_MESSAGE(about, /INF, TITLE='About FEATURETOOL')
   end
   
   'QUIT': begin
      WIDGET_CONTROL, ev.TOP, /DESTROY        
      return
   end
   
   ELSE: print, 'Event type not yet implemented'
   
ENDCASE

; Process the present file according to the present settings
if s.havefile and doimage then begin
   if doclearhist then s.p.hist *= 0L
   f = featuretool_processimage(s)
   if savefeatures then write_gdf, f, fn
endif

; return the present state of all variables
WIDGET_CONTROL, ev.TOP, SET_UVALUE=s

end

;;;
;;; FEATURETOOL
;;;
;;; The main routine
;;;
pro featuretool, a

COMPILE_OPT IDL2

w = 640                         ; default image dimensions
h = 480
file = ''                       ; current file name
havefile = 0
if n_params() eq 1 then begin
   sz = size(a, /dimensions)
   w = sz[0]
   h = sz[1]
   file = filepath("featuretool.pgm",/TMP)
   havefile = 1
   write_ppm, file, a
endif

hist = lonarr(50)                ; histogram of sub-pixel locations

; p contains all the parameters used or reported by BPASS and FEATURE.
; Storing them here avoids the need for common blocks, and thus allows
; for multiple instances of FEATURETOOL. 
p = { xmin:0, $
      xmax:w-1, $
      ymin:0, $
      ymax:h-1, $
      w:w, $
      h:h, $
      deinterlace:0, $
      invert:0, $
      noise:1., $
      extent:7, $
      diameter:7, $
      separation:8, $
      min:1, $
      masscut:0, $
      nlmax:0, $
      nfeatures:0, $
      hist:hist}

; The base widget: contains the controls and the tabbed windows.
base = WIDGET_BASE(TITLE = 'Feature Tool', MBAR = bar, UVALUE = 'base', /COLUMN)

file_menu = WIDGET_BUTTON(bar, VALUE = 'File', /MENU)
file_bn = WIDGET_BUTTON(file_menu, VALUE = 'Save &Image...', $
                        UVALUE = 'SAVEIMAGE',  ACCELERATOR = 'Alt+I')
file_bn = WIDGET_BUTTON(file_menu, VALUE = '&Save Features...', $
                         UVALUE = 'SAVEFEATURES', ACCELERATOR = 'Alt+S')
file_bn = WIDGET_BUTTON(file_menu, VALUE = 'Save &Routine...', $
                         UVALUE = 'SAVEPRO', ACCELERATOR = 'Alt+R')

sensitive = strlen(file_which("idlsnap.pro")) gt 0
file_bn = WIDGET_BUTTON(file_menu, VALUE = 'IDLSNA&P', $
                        UVALUE = 'IDLSNAP', ACCELERATOR = 'Alt+P', $
                        SENSITIVE = sensitive)

file_bn = WIDGET_BUTTON(file_menu, VALUE = '&Quit', $
                        UVALUE = 'QUIT', ACCELERATOR = 'Alt+Q')

help_menu = WIDGET_BUTTON(bar, VALUE = 'Help', /MENU, /ALIGN_RIGHT)
help_bn1 = WIDGET_BUTTON(help_menu, VALUE = 'Instructions', UVALUE = 'HELPINSTRUCTIONS')
help_bn2 = WIDGET_BUTTON(help_menu, VALUE = 'About', UVALUE = 'HELPABOUT')

; All controls (file selection, preprocessing, bpass, feature, and
; reporting) are contained in this widget.
controlbar = WIDGET_BASE(base, column=2)

; File selector for image files.
fn =  CW_FILESEL(controlbar, FILTER = ['.pgm, .ppm', 'png', 'jpg'], $
                 /IMAGE_FILTER, UVALUE = 'FILE')

; preprossing, bpass, feature and reporting widgets
inputsbar = WIDGET_BASE(controlbar, column = 1)

; First line: region of interest
roibar = WIDGET_BASE(inputsbar, column = 4)
wxmin = CW_FIELD(roibar, title='XMIN', UVALUE='XMIN', $
                 value = p.xmin, /INTEGER, /RETURN_EVENTS, /COLUMN)
wxmax = CW_FIELD(roibar, title='XMAX', UVALUE='XMAX', $
                 value = p.xmax, /INTEGER, /RETURN_EVENTS, /COLUMN)
wymin = CW_FIELD(roibar, title='YMIN', UVALUE='YMIN', $
                 value = p.ymin, /INTEGER, /RETURN_EVENTS, /COLUMN)
wymax = CW_FIELD(roibar, title='YMAX', UVALUE='YMAX', $
                 value = p.ymax, /INTEGER, /RETURN_EVENTS, /COLUMN)

; Second line: preprossing and bpass
bpassbar = WIDGET_BASE(inputsbar, column = 3, /ALIGN_RIGHT)
; preprocessing
wbppreprocess = CW_BGROUP(bpassbar, ['Deinterlace', 'Invert'], $
                          BUTTON_UVALUE = ['BPDEINTERLACE', 'BPINVERT'], $
                          COLUMN = 2, /NONEXCLUSIVE, UVALUE = 'BPPREPROCESS', $
                          LABEL_TOP = 'Preprocess', /FRAME, XSIZE = 275)
; bpass: noise
wbpnoise = CW_FIELD(bpassbar, title='NOISE', $
                    UVALUE='BPNOISE', value = p.noise, $
                    /FLOATING, /RETURN_EVENTS, /COLUMN)
; bpass: extent
wbpextent = CW_FIELD(bpassbar, title='EXTENT', $
                     UVALUE='BPEXTENT', value = p.extent, $
                     /INTEGER, /RETURN_EVENTS, /COLUMN)

; Third line: feature
featurebar = WIDGET_BASE(inputsbar, column = 4)
; feature: particle diameter
wfdiameter = CW_FIELD(featurebar, title='DIAMETER', $
                      UVALUE='FDIAMETER', value = p.diameter, $
                      /INTEGER, /RETURN_EVENTS, /COLUMN)
; feature: particle separation
wfseparation = CW_FIELD(featurebar, title='SEPARATION', $
                        UVALUE='FSEPARATION', value = p.separation, $
                        /INTEGER, /RETURN_EVENTS, /COLUMN)
; feature: minimum brightness to consider
wfmin = CW_FIELD(featurebar, title='MIN', $
                 UVALUE='FMIN', value = p.min, $
                 /INTEGER, /RETURN_EVENTS, /COLUMN)
; feature; minimum integrated brightness to retain
wfmasscut = CW_FIELD(featurebar, title='MASSCUT', $
                     UVALUE='FMASSCUT', value = p.masscut, $
                     /INTEGER, /RETURN_EVENTS, /COLUMN)

; Fourth line: Results
resultsbar = WIDGET_BASE(inputsbar, column = 2, TITLE = 'Results', $
                         /ALIGN_CENTER, /FRAME)
; results: number of local maxima
wlmax = CW_FIELD(resultsbar, title = 'Local Maxima', value = p.nlmax, $
                 /INTEGER, /NOEDIT, /COLUMN)
; results: number of features
wfeat = CW_FIELD(resultsbar, title = 'Features', value = p.nfeatures, $
                 /INTEGER, /NOEDIT, /COLUMN)

; Tabbed windows showing results
wtab = WIDGET_TAB(base)
; plot of image with overlaid features
wimagetab = WIDGET_BASE(wtab, TITLE = 'IMAGE', /COLUMN)
wimage = WIDGET_DRAW(wimagetab, xsize = 640, ysize = 480, RETAIN = 2)
; plot of image filtered by BPASS
wbpasstab = WIDGET_BASE(wtab, TITLE = 'BPASS', /COLUMN)
wbpass = WIDGET_DRAW(wbpasstab, xsize = 640, ysize = 480, RETAIN = 2)
; plot of features' radius of gyration versus integrated brightness
wmasscuttab = WIDGET_BASE(wtab, TITLE = 'MASSCUT', /COLUMN)
wmasscut = WIDGET_DRAW(wmasscuttab, xsize = 640, ysize = 480, RETAIN = 2)
; plot of fractional feature locations (pixel bias detection)
wdiametertab = WIDGET_BASE(wtab, TITLE = 'DIAMETER', /COLUMN)
wdiameter = WIDGET_DRAW(wdiametertab, xsize = 640, ysize = 480, RETAIN = 2)

; structure containing identifiers for all the widgets
w = {wxmin:wxmin, wxmax:wxmax, wymin:wymin, wymax:wymax, $
     wimage:wimage, wbpass:wbpass, wmasscut:wmasscut, wdiameter:wdiameter, $
     wtab:wtab, wlmax:wlmax, wfeat:wfeat}

; structure containing the widget identifiers, information on the
; current file and parameters used and reported by BPASS and FEATURE.
s = {w:w, file:file, havefile:havefile, p:p}

WIDGET_CONTROL, base, /REALIZE

if havefile then $
   f = featuretool_processimage(s)

WIDGET_CONTROL, base, SET_UVALUE = s
XMANAGER, 'featuretool', base, /NO_BLOCK
end
