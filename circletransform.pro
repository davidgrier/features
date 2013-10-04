;+
; NAME:
;    circletransform
;
; PURPOSE:
;    Performs a transform similar to a Hough transform
;    for detecting circular features in an image.
;
; CATEGORY:
;    Image analysis, feature detection
;
; CALLING SEQUENCE:
;    b = circletransform(a)
;
; INPUTS:
;    a: [nx,ny] image data
;
; KEYWORD PARAMETERS:
;    noise: estimate for additive pixel noise.
;        Default: noise estimated by MAD().
;    deinterlace: if set to an odd number, then only perform
;        transform on odd field of an interlaced image.
;        If set to an even number, transform even field.
;        Default: Not set or set to zero: transform entire frame.
;
; OUTPUTS:
;    b: [nx,ny] circle transform.  Peaks correspond to estimated
;        centers of circular features in a.
;
; KEYWORD OUTPUTS:
;    range: mean range used in tallying votes.
;
; PROCEDURE:
;    Compute the gradient of the image.  The local gradient at each
;    pixel defines a line along which the center of a circle may
;    lie.  Cast votes for pixels along the line in the transformed
;    image.  The pixels in the transformed image with the most votes
;    correspond to the centers of circular features in the original
;    image.
;
; REFERENCE:
; F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao, L. Dixon
; & D. G. Grier,
; Flow visualization and flow cytometry with holographic video
; microscopy, Optics Express 17, 13071-13079 (2009)
;
; EXAMPLE:
;    IDL> b = circletransform(a)
;
; MODIFICATION HISTORY:
; 10/07/2008 Written by David G. Grier, New York University.
; 01/26/2009 DGG Added DEINTERLACE keyword. Gracefully handle
;    case when original image has no features. Documentation cleanups.
; 02/03/2009 DGG Replaced THRESHOLD keyword with NOISE.
; 06/10/2010 DGG Documentation fixes.  Added COMPILE_OPT.
; 05/02/2012 DGG Updated keyword parsing.  Formatting.
; 06/24/2012 DGG Streamlined index range checking in inner loop
;    to improve efficiency.
; 07/16/2012 DGG IMPORTANT: Center results on pixels, not on vertices!
;    Use array_indices for clarity.
; 11/10/2012 DGG Default range should be an integer.
;    Returned array should be cast to integer, not float
; 11/23/2012 DGG Use Savitzky-Golay estimate for derivative.
;    Eliminate SMOOTHFACTOR parameter.  Upgrade parameter checking.
; 11/25/2012 DGG Limit search range by uncertainty in gradient
;    direction.  Remove RANGE keyword.
; 11/30/2012 DGG Optionally return mean range as RANGE keyword.
; 01/16/2013 DGG estimate noise with MAD() by default.
; 01/24/2013 DGG correct test for deinterlace = 0.
; 02/09/2013 DGG use savgol2d() to compute derivatives.
;    Displace by half a pixel to center, not a whole pixel.
; 02/17/2013 DGG RANGE is the median range of voting pixels, not the
;    mean.
; 03/04/2013 DGG shift by +1 rather than by +0.5.  Limit range if
;    noise is very small.
; 03/17/2013 DGG calculate coordinates explicitly rather than using
;    array_indices, which turns out to be slow.  More efficient array
;    indexing.  No more need to shift pixels for alignment.
; 03/27/2013 DGG eliminate repetitive operations in loops.
; 05/13/2013 DGG suppress borders, which are over-counted.
; 10/04/2013 DGG and Mark Hannel: fix boundary cropping.
;
; Copyright (c) 2008-2013 David G. Grier and Mark Hannel
;
;-

function circletransform, a_, $
                          noise = noise, $
                          range = range, $
                          deinterlace = deinterlace

COMPILE_OPT IDL2

umsg = 'USAGE: b = circletransform(a)'

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

dodeinterlace = isa(deinterlace, /scalar, /number) ? deinterlace gt 0 : 0
if dodeinterlace then begin
   n0 = deinterlace mod 2
   a = float(a_[*, n0:*:2])
endif else $
   a = float(a_)

if ~isa(noise, /scalar, /number) then $
   noise = mad(a)

dx = savgol2d(7, 3, dx = 1)
dadx = convol(a, dx, /edge_truncate)
dady = convol(a, transpose(dx), /edge_truncate)
if dodeinterlace then dady /= 2.
grada = sqrt(dadx^2 + dady^2)           ; magnitude of the gradient
dgrada = noise * sqrt(2. * total(dx^2)) ; error in gradient magnitude due to noise
w = where(grada gt 2.*dgrada, npts)     ; only consider votes with small angular uncertainty

b = intarr(nx, ny)              ; accumulator array for the result

if npts le 0 then return, b

xp = w mod nx                   ; coordinates of pixels with strong gradients
yp = w / nx
if dodeinterlace then yp = 2*yp + n0

grada = grada[w]                ; gradient direction at each pixel
dgrada = dgrada[w] / grada
costheta = dadx[w] / grada
sintheta = dady[w] / grada

rng = round(2./tan(dgrada/2.)) < nx
range = max(rng)
r = findgen(2*range + 1) - range

nx--
ny--
for i = 0L, npts-1L do begin 
   rr = r[range-rng[i]:range+rng[i]]
   x = (xp[i] + round(rr * costheta[i])) > 0 < nx
   y = (yp[i] + round(rr * sintheta[i])) > 0 < ny
   b[x, y]++
endfor

; borders are over-counted because of > and <
b[*, 0] = 0
b[0, *] = 0
b[-1, *] = 0
b[*, -1] = 0

if arg_present(range) then $
   range = median(rng)

return, b
end
