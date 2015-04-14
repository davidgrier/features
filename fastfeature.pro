;+
; NAME:
;    fastfeature
;
; PURPOSE:
;    Find the centroids of disk-like features in an image quickly.
;
; CATEGORY:
;    Image processing; digital video microscopy
;
; CALLING SEQUENCE:
;    f = fastfeature(image,threshold)
;
; INPUTS:
;    image: two-dimensional gray-scale image
;
;    threshold: gray level distinguishing objects from background.
;
; KEYWORD FLAGS:
;    dark: if set, then objects are darker than background.
;        Default: objects are assumed to be brighter than
;        background.
;
;    npixels: if set, return number of pixels in each region
;        Default: return integrated brightness relative to threshold.
;
;    center: if set, return coordinates relative to the center of
;        the field of view.
;        Default: origin is set at the lower left corner.
;
; KEYWORD PARAMETERS:
;    deinterlace: if set to an even or odd integer then analyze only
;        the even or odd field in an interlaced image, respectively.
;        Default: not set, or set to 0: analyze entire image.
;
;    pickn: return the brightest N features
;
;    count: number of features returned.
;
; OUTPUTS:
;    f: [nd+1,nobjects] array of located features.
;        f[0:nd-1,*]: position in nd dimensions
;        f[nd,*]: integrated brightness of the feature relative to threshold
;
; RESTRICTIONS:
;    Images should be two dimensional with real-valued pixels.
;
; PROCEDURE:
;    Threshold the image.  Locate contiguous blobs of
;    above-threshold pixels.  Return the brightness-weighted
;    centroids of those blobs.
;
; EXAMPLE:
;    IDL> a = idlsnap()   ; acquire image
;    IDL> h = histogram(a,min=0,max=255)
;    IDL> plot, h         ; estimate threshold from plot
;    IDL> f = fastfeature(a,threshold)
;    IDL> foverlay, f, a  ; see how it turned out
;
; MODIFICATION HISTORY:
; 12/12/2004 David G. Grier, New York University: created.
; 01/14/2009 DGG return -1 if no particles are found above threshold.
; 06/10/2010 DGG Documentation fixes.  Added COMPILE_OPT.
; 05/12/2012 DGG Added PICKN for compatibility with feature.
; 06/01/2012 David B. Ruffner, NYU: Generalized to nd dimensions
;    for images (nd = 2), volumetric data (nd = 3), and more
;    general data sets.
; 06/22/2012 DGG Improved input checking.
; 07/18/2012 DGG Fixed half-pixel offset with CENTER keyword.
;    Added DEINTERLACE keyword.
; 10/16/2012 DGG Added COUNT keyword.
; 10/25/2012 DGG COUNT should not include background as a feature
;    (duh).
; 02/17/2013 DGG Setting DEINTERLACE = 0 does not deinterlace.
; 02/26/2013 DGG Clean up threshold code.  Set ALL_NEIGHBORS for
;    label_region.
; 03/17/2013 DGG More efficient array indexing.  Simplified
;    deinterlace code.  Simplified main loop.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 05/12/2013 DGG fix brightness weighting.  Added NPIXELS.
;
; Copyright (c) 2004-2013 David G. Grier and David B. Ruffner
;-

function fastfeature, image, threshold, $
                      center = center, $
                      dark = dark, $
                      pickn = pickn, $
                      count = count, $
                      deinterlace = deinterlace, $
                      npixels = npixels

  COMPILE_OPT IDL2

  umsg = 'p = fastfeature(image, threshold)'

  if n_params() ne 2 then begin
     message, umsg, /inf
     return, -1
  endif

  if ~isa(image, /number, /array) then begin
     message, umsg, /inf
     message, 'image must be a numerical array', /inf
     return, -1
  endif

  sz = size(image)
  nd = sz[0]

  if nd ne 2 && nd ne 3 then begin
     message, umsg, /inf
     message, 'image must be two- or three-dimensional', /inf
     return, -1
  endif

  if ~isa(threshold, /number, /scalar) then begin
     message, umsg, /inf
     message, 'threshold must be a number', /inf
     return, -1
  endif
   
  dodeinterlace = keyword_set(deinterlace)
  countpixels = keyword_set(npixels)
  y0 = (dodeinterlace) ? long(deinterlace) mod 1L : 0
  dy = (dodeinterlace) ? 2 : 1
  img = keyword_set(dark) ? threshold - image[*, y0:*:dy, *] : $
        image[*, y0:*:dy, *] - threshold
  reg = label_region(img gt 0, /all_neighbors)

  ;;; Find centroid of each labeled region
  n = histogram(reg, reverse_indices = r)
  count = n_elements(n) - 1     ; do not count background as a feature
  if count le 0 then $
     return, -1

  f = fltarr(nd+1, count, /nozero)
  for i = 1, count do begin       ; the background is region 0
     ndx = r[r[i]:r[i+1]-1]       ; 1D indices of pixels in region i
     nn = array_indices(img, ndx) ; nd-dimensional indices of pixels in i
     if dodeinterlace then $
        nn[1,*] = dy*temporary(nn[1,*]) + y0
     v = transpose(rebin(img[ndx], n[i], nd, /sample))        ; values in region i
     f[0,i-1] = (n[i] eq 1) ? nn : total(nn*v, 2)/total(v, 2) ; value-weighted centers
     f[nd,i-1] = (countpixels) ? n[i] : total(img[ndx])       ; integrated brightness
  endfor

  if isa(pickn, /scalar, /number) then begin
     if count gt pickn then begin
        order = reverse(sort(f[nd,*]))
        f = f[*, order[0:pickn-1]]
        count = pickn
     endif
  endif

  if keyword_set(center) then begin
     for i = 0, nd-1 do $
        f[i, *] -= (sz[i+1] - 1.)/2.
  endif
  
  return, f
end
