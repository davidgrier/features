;+
; NAME:
;    bpass
;
; PURPOSE:
;    Implements a real-space bandpass filter which suppress 
;    pixel noise and long-wavelength image variations while 
;    retaining information about features of a characteristic size.
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    res = bpass(image, lshort, llong)
; 
; INPUTS:
;    image: The two-dimensional array to be filtered.
;    lshort: Short-wavelength cut-off.
;    llong: Long-wavelength cut-off.
;
; OUTPUTS:
;    res: filtered image.
;
; PROCEDURE:
;    Simple 'wavelet' convolution.
;
; MODIFICATION HISTORY:
; Written by David G. Grier, The University of Chicago, 2/93.
; 05/95 DGG Greatly revised version: 1D convolutions rather than 2D
; 12/95 John C. Crocker.  Added /field keyword.
; 08/99 DGG Memory optimizations and fixed normalization.
; 09/01/1999 DGG. Added EDGE_TRUNCATE to save a little more data.
; 09/10/2010 David G. Grier, New York University.  Added COMPILE_OPT
;    and revised documentation.
;
; Copyright 1993-2010 David G. Grier and John C. Crocker
;-
function bpass, image, lshort, llong, $
                field = field, noclip = noclip
COMPILE_OPT IDL2
;on_error, 2				; go to caller on error

b = float(lshort)
w = round(llong > (2. * b))
N = 2*w + 1

; Gaussian convolution kernel: removes short-wavelength noise
r = (findgen(N) - w)/(2. * b)
gx = exp(-r^2) / (2. * b * sqrt(!pi))
gy = transpose(gx)

; Boxcar average kernel: calculates background
bx = fltarr(N) + 1./N
by = transpose(bx)

if keyword_set(field) then begin
   if N mod 4 eq 1 then $
      ndx = 2*indgen(w+1)$
   else $
      ndx = 1 + 2*indgen(w)
   gy = 2. * gy(ndx)           ; fixes normalization
   nn =  n_elements(ndx)
   by = fltarr(nn) + 1./nn
endif

res = float(image)
; Gaussian convolutions
g = convol(res, gx, /edge_truncate)
g = convol(temporary(g), gy, /edge_truncate)
; Smoothed background.  Note that we can't use SMOOTH
; in case analysis is field-based, rather than frame-based.
res = convol(temporary(res), bx, /edge_truncate)
res = convol(temporary(res), by, /edge_truncate)
res = g - temporary(res)

if keyword_set(noclip) then $
   return, res $
else $
   return, res > 0

end

;*********** end of bpass.pro












