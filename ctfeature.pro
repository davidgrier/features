;+
; NAME:
;    ctfeature
;
; PURPOSE:
;    Identify ring-like features in images
;
; CATEGORY:
;    Image analysis, feature detection
;
; CALLING SEQUENCE:
;    f = ctfeature(a)
;
; INPUTS:
;    a: two-dimensional image data
;
; KEYWORD PARAMETERS:
;    noise: estimate for additive noise in pixel values
;        Default: set by CIRCLETRANSFORM
;
;    range: range over which to search for features [pixels]
;        Default: set by CIRCLETRANSFORM.
;
;    threshold: threshold for detecting features
;        Default: estimated from range computed by CIRCLETRANSFORM. 
;
;    pickn: number of features to seek, brightest first
;        Default: all
;
;    count: number of features returned.
;
; KEYWORD FLAGS:
;    deinterlace: Set to an even number to find features
;        in the even field, odd in the odd.
;
;    quiet: If set, do not print informational messages.
;
; OUTPUTS:
;    f: [2,npts] array of feature coordinates
;
; PROCEDURE:
;    CIRCLETRANSFORM transforms ring-like features in an image into
;        bright features on a dark background.
;    FASTFEATURE locates these bright features.
;
; MODIFICATION HISTORY:
; 10/15/2012 Written by David G. Grier, New York University
; 11/10/2012 DGG Added QUIET keyword.  Changed default smooth factor
;   from 3 to 5.  Cast threshold to an integer.
; 11/23/2012 DGG Updated for consistency with CIRCLETRANSFORM.
;   Removed SMOOTHFACTOR parameter.
; 11/25/2012 DGG Removed NOISE and RANGE parameters.
; 12/21/2012 DGG Pass NOISE estimate to CIRCLETRANSFORM to take
;   advantage of new range-estimation code.
; 01/16/2013 DGG accept RANGE keyword.
;
; Copyright (c) 2012 David G. Grier
;
;-
function ctfeature, a, $
                    noise = noise, $
                    threshold = threshold, $
                    pickn = pickn, $
                    count = count, $
                    ct = ct, $
                    deinterlace = deinterlace, $
                    quiet = quiet

COMPILE_OPT IDL2

umsg = 'USAGE: f = ctfeature(a)'

if n_params() ne 1 then begin
   message, umsg, /inf
   return, -1
endif

noprint = keyword_set(quiet)

; Find candidate features ...
;; transform ring-like patterns into spots
ct = circletransform(a, noise = noise, range = range, deinterlace = deinterlace)

;; centers of spots are estimates for particle centers: (xp, yp)
if ~isa(threshold, /number, /scalar) then begin
   threshold = round(!pi * range^2 / 4.)
   if keyword_set(deinterlace) then threshold /= 2
endif

p = fastfeature(ct, threshold, pickn = pickn, count = count) ; find peaks

if count lt 1 then begin
   message, umsg, /inf, noprint = noprint
   message, 'no features found above threshold = '+strtrim(threshold, 2), /inf, noprint = noprint
   return, -1
endif

return, p
end
