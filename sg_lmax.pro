;+
; NAME:
;    sg_lmax
;
; PURPOSE:
;    Find local maxima in two-dimensional data sets using
;    one-dimensional Savitzky-Golay filters.
;
; CATEGORY:
;    Image processing
;
; CALLING SEQUENCE:
;    r = sg_lmax(image, extent, [separation], 
;
; INPUTS:
;    image: two-dimensional numerical data set
;    extent: diameter of features
;
; OPTIONAL INPUTS:
;    separation: minimum allowable separation between
;        features' centers.  Default value is extent+1
;
; KEYWORD PARAMETERS:
;    min: Minimum value to be considered for a local maximum.
;        Default: min = 0.
;
; KEYWORD FLAGS:
;    field: If set, scale by 2 in the y direction.
;
; OUTPUTS:
;    r: [2,npts] array of (x,y) coordinates of local maxima.
;
; MODIFICATION HISTORY:
; 07/29/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-
function sg_lmax, a, extent, separation, min = min, field = field

COMPILE_OPT IDL2

min = isa(min, /scalar, /number) ? float(min) : 0.

order = 5
sx = savgol(extent, extent, 0, order)
kx = savgol(separation, separation, 1, order)
if keyword_set(field) then begin
   sy = transpose(savgol(extent/2, extent/2, 0, order))
   ky = transpose(savgol(separation/2, separation/2,  1,  order))
endif else begin
   sy = transpose(sx)
   ky = transpose(kx)
endelse
b = float(a)
bx = convol(convol(b, sy, /edge_truncate), kx, /edge_truncate)
by = convol(convol(b, sx, /edge_truncate), ky, /edge_truncate)

c = (bx*shift(bx, 1, 0) lt 0) and (by*shift(by, 0, 1) lt 0)
w = where(c and (b gt min), count)
r = array_indices(a, w)
return, r
end
