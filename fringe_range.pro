;+
; NAME:
;    fringe_range
;
; PURPOSE:
;    Identify range of pixels in an image of circularly symmetric
;    oscillatory features such as holographic fringe patterns.
;
; CATEGORY:
;    Image analysis, feature detection
;
; CALLING SEQUENCE:
;    b = fringe_range(a, p)
;
; INPUTS:
;    a: [nx,ny] gray-scale image data
;    p: [2, npts] coordinates of target centers
;
; KEYWORD PARAMETERS:
;    nfringes: Number of extrema to retain.
;        Default: 20 (10 bright rings).
;
;    deinterlace: if set to an odd number, then only perform
;        transform on odd field of an interlaced image.
;        not set or 0: analyze entire frame
;        odd: analyze odd field
;        even: analyze even field
;
; OUTPUTS:
;    rad: [npts] array of range estimates for each input point.
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
; 08/05/2013 Written by David G. Grier, New York University.
; 08/06/2013 Skip cases with no extrema.
;
; Copyright (c) 2013 David G. Grier
;-

function fringe_range, a, rp, $
                       nfringes = nfringes, $
                       deinterlace = deinterlace

COMPILE_OPT IDL2

umsg = 'USAGE: rad = fringe_range(image, centers)'

if n_params() ne 2 then begin
   message, umsg, /inf
   return, -1
endif

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   return, -1
endif
if size(a, /n_dimensions) ne 2 then begin
   message, umsg, /inf
   message, 'IMAGE must be a two-dimensional numeric array', /inf
   return, -1
endif
sz = size(a, /dimensions)
nx = sz[0]
ny = sz[1]

if ~isa(rp, /number, /array) then begin
   message, umsg, /inf
   message, 'CENTERS must be an array of target coordinates', /inf
   return, -1
endif
sz = size(rp)
if sz[0] gt 2 or sz[1] lt 2 then begin
   message, umsg, /inf
   message, 'CENTERS should be a [2, nfeatures] array of coordinates', /inf
   return, -1
endif
nfeatures = (sz[0] eq 2) ? sz[2] : 1

nfringes = (isa(nfringes, /scalar, /number)) ? long(nfringes) : 20

rad = intarr(nfeatures)         ; the answer

for ndx = 0, nfeatures - 1 do begin
   rc = rp[0:1, ndx]
   aa = aziavg(a, center = rc, deinterlace = deinterlace)
   rn = extrema(aa, count = count)
   if count le 0 then continue
   n = nfringes < (count - 1)
   rad[ndx] = rn[n]
endfor

return, rad
end
