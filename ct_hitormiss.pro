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
;    p: (x, y)  target location to hit
;
; KEYWORD PARAMETERS:
;    range: range over which a circle's center will be sought.
;        Default: 100 pixels
;
;    noise: estimate for additive pixel noise.
;        Default: noise estimated by MAD().
;
;    deinterlace: if set to an odd number, then only perform
;        transform on odd field of an interlaced image.
;        If set to an even number, transform even field.
;
; KEYWORD FLAGS:
;    coordinates: If set, return coordinates of hit/miss pixels
;        [x,y,hm] where hm is 1 for hits and 0 for misses
;
;    distance: If set, return image with contributing pixels labeled
;        by their distance to the target point.
;        If COORDINATES also is set, then return
;        [x,y,dist]
;
; OUTPUTS:
;    b: [nx,ny] circle transform.  Peaks correspond to estimated
;        centers of circular features in a.
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
; 01/16/2012 DGG estimate noise with MAD() by default.
;
; Copyright (c) 2012-2013 David G. Grier
;-

function ct_hitormiss, a_, p, $
                       range = range, $
                       noise = noise, $
                       deinterlace = deinterlace, $
                       distance = distance, $
                       coordinates = coordinates

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

if ~isa(range, /scalar, /number) then range = 100

dodeinterlace = isa(deinterlace, /scalar, /number)
if dodeinterlace then begin
   n0 = deinterlace mod 2
   a = float(a_[*, n0:*:2])
endif else $
   a = float(a_)

if ~isa(noise, /scalar, /number) then $
   noise = mad(a)

returncoordinates = keyword_set(coordinates)

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
w = where(grada gt 2.*dgrada, npts)     ; select points with small angular uncertainty

hit = 0.*a

if npts le 0 then $
   return, (returncoordinates) ? -1 : hit

xy = array_indices(a, w)
if dodeinterlace then xy[1,*] = 2.*xy[1,*] + n0
xy += 1.

qx = xy[0,*] - p[0]
qy = xy[1,*] - p[1]

ww = where((qx^2 + qy^2) le range^2, npts)
if npts le 0 then $
   return, (returncoordinates) ? -1 : hit

w = w[ww]
qx = qx[ww]
qy = qy[ww]
grada = grada[w]
costheta = dadx[w] / grada
sintheta = dady[w] / grada
delta = abs(qx * sintheta - qy * costheta)
ddelta = dgrada * abs((qx * costheta + qy * sintheta)) / grada

if keyword_set(distance) then begin
   if returncoordinates then $
      return, [xy[*,ww], transpose(delta)]
   hit[w] = delta
   return, hit
endif

if returncoordinates then $
   return, [xy[*, ww], transpose(delta) le ddelta]

hit[w] = -1.                    ; assume all points have missed

ww = where(delta le ddelta, npts)
if npts le 0 then return, hit

hit[w[ww]] = 1.
return, hit

end
