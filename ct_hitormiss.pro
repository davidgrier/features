;+
; NAME:
;    ct_hitormiss
;
; PURPOSE:
;    Identify pixels in an image that are transformed to a particular
;    feature by circletransform
;
; CATEGORY:
;    Image analysis, feature detection
;
; CALLING SEQUENCE:
;    b = ct_hitormiss(a, p)
;
; INPUTS:
;    a: [nx,ny] gray-scale image data
;    p: [2,npts]  target locations to check for hits
;
; KEYWORD PARAMETERS:
;    noise: estimate for additive pixel noise.
;        Default: noise estimated by MAD().
;
;    deinterlace: if set to an odd number, then only perform
;        transform on odd field of an interlaced image.
;        If set to an even number, transform even field.
;
; OUTPUTS:
;    b: [nx,ny] Hit or miss map
;       values:
;       -1: miss for all features
;       0: no determination
;       n: hit for feature n
;
; PROCEDURE:
;    Compute the gradient of the image.  The local gradient at each
;    pixel defines a line along which the center of a circle may
;    lie.  Pixel "hits" if the line passes within range of the
;    specified feature.
;
; REFERENCE:
; F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao, L. Dixon
; & D. G. Grier,
; Flow visualization and flow cytometry with holographic video
; microscopy, Optics Express 17, 13071-13079 (2009)
;
; MODIFICATION HISTORY:
; 05/02/2012 Written by David G. Grier, New York University.
; 07/16/2012 DGG Use array_indices for clarity.
; 07/18/2012 DGG coordinates return 0 for miss rather than -1.
; 11/23/2012 DGG Update for consistency with CIRCLETRANSFORM.
;    Gradients estimated with Savitzky-Golay convolution kernel.
;    Eliminated SMOOTHFACTOR parameter.  Precision set by estimated
;    error in voting direction.  Eliminated PRECISION keyword.
; 11/30/2012 DGG Fixed hit test for /coordinates.
; 01/16/2013 DGG estimate noise with MAD() by default.
; 01/21/2013 DGG complete overhaul: hit or miss for multiple target
;    points.  Remove COORDINATES and DISTANCE keywords.
; 01/22/2013 DGG Use CLUSTER() to restrict hits nearest-neighborhoods.
;    Removed RANGE keyword.
;
; Copyright (c) 2012-2013 David G. Grier
;-

function ct_hitormiss, a_, p, $
                       noise = noise, $
                       deinterlace = deinterlace

COMPILE_OPT IDL2

umsg = 'USAGE: b = ct_hitormiss(a, p)'

if ~isa(a_, /number, /array) then begin
   message, umsg, /inf
   return, -1
endif
if size(a_, /n_dimensions) ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional numeric array', /inf
   return, -1
endif
sz = size(a_, /dimensions)
nx = sz[0]
ny = sz[1]

if ~isa(p, /number, /array) then begin
   message, umsg, /inf
   message, 'P must be an array of target coordinates', /inf
   return, -1
endif
sz = size(p)
if sz[0] gt 2 or sz[1] lt 2 then begin
   message, umsg, /inf
   message, 'P should be a [2,npts] array of coordinates', /inf
   return, -1
endif
npts = (sz[0] eq 2) ? sz[2] : 1

dodeinterlace = isa(deinterlace, /scalar, /number)
if dodeinterlace then begin
   n0 = deinterlace mod 2
   a = float(a_[*, n0:*:2])
endif else $
   a = float(a_)

if ~isa(noise, /scalar, /number) then $
   noise = mad(a)

hit = 0*a

; Third-order two-dimensional Savitzky-Golay filter over 5x5 image patch
dx = [[ 0.0738, -0.1048,  0.0000,  0.1048, -0.0738], $
      [-0.0119, -0.1476,  0.0000,  0.1476,  0.0119], $
      [-0.0405, -0.1619,  0.0000,  0.1619,  0.0405], $
      [-0.0119, -0.1476,  0.0000,  0.1476,  0.0119], $
      [ 0.0738, -0.1048,  0.0000,  0.1048, -0.0738]]
dadx = convol(a, dx, /center, /edge_truncate)
dady = convol(a, transpose(dx), /center, /edge_truncate)
if dodeinterlace then dady /= 2.
grada = sqrt(dadx^2 + dady^2)           ; magnitude of the gradient
dgrada = noise * sqrt(2. * total(dx^2)) ; error in gradient estimate due to noise

w = where(grada gt 2.*dgrada, ngood)    ; select points with small angular uncertainty
if ngood le 0 then $
   return, hit

grada = grada[w]
dgrada = dgrada[w]/grada
costheta = dadx[w]/grada
sintheta = dady[w]/grada

hit[w] = -1                     ; all points start out as misses

xy = array_indices(a, w)
if dodeinterlace then xy[1,*] = 2.*xy[1,*] + n0
xy += 1. ; is this needed for the /center flag on dadx, dady?

id = (npts gt 1) ? cluster(xy, p[0:1, *]) : intarr(npts)

for n = 0, npts-1 do begin
   qx = xy[0,*] - p[0, n]
   qy = xy[1,*] - p[1, n]
   rsq = qx^2 + qy^2
   delta = abs(qx * sintheta - qy * costheta)
   ddelta = dgrada * abs((qx * costheta + qy * sintheta)) < 5.
   ww = where((delta le ddelta) and (id eq n), nhit)
   if nhit gt 0 then $
      hit[w[ww]] = n + 1
endfor

return, hit

end
