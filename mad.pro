;+
; NAME:
;    mad
;
; PURPOSE:
;    Computes median absolute difference (MAD) estimate for the noise
;    in an image.
;
; CATEGORY:
;    Image processing
;
; CALLING SEQUENCE:
;    noise = mad(a)
;
; INPUTS:
;    a: image data
;
; KEYWORD PARAMETERS:
;    w: width of median window.  Default: 10 [pixels]
;
; OUTPUTS:
;    noise: estimate for the amplitude of the Gaussian random noise
;           in the image.
;
; MODIFICATION HISTORY:
; 01/15/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-
function mad, a, w = w

COMPILE_OPT IDL2

umsg = 'noise = mad(a)'

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   return, -1.
endif

if ~isa(w, /number, /scalar) then $
   w = 10

return, median(abs(float(a) - median(a, w)))
end
