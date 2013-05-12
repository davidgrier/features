;+
; NAME:
;    Feature	
;
; PURPOSE:
;    Finds and measures roughly circular 'features' within an image.
;
; CATEGORY:
;    Image Analysis
;
; CALLING SEQUENCE:
;    f = feature(image, diameter, [separation])
;
; INPUTS:
;    image: (nx,ny) array that presumably contains some
;        features worth finding
;    diameter: a parameter that should be a little greater than
;        the diameter of the largest features in the image.
;        NOTE: Diameter MUST BE ODD valued.
;
; OPTIONAL INPUTS:
;    separation: an optional parameter that specifies the 
;        minimum allowable separation between feature 
;        centers. The default value is diameter+1.
;
; KEYWORD PARAMETERS:
;    min: Minimum pixel intensity considered to be a
;        candidate local maximum (default: 64th percentile).
;        Setting min appropriately improves performance.
;    masscut: Discard candidate features with integrated
;        brightness below this value.
;        Setting masscut appropriately improves performance.
;    pickn: Return only the PICKN brightest candidate features.
;    maxits: Maximum number of iterations (default = 10).
;
; KEYWORD FLAGS:
;    field:  Set this keyword if image is actually just one field
;        of an interlaced (e.g. video) image. All the masks
;        will then be constructed with a 2:1 aspect ratio.
;    iterate: if the refined centroid position is too far from
;        the initial estimate, iteratively recalculate the centroid
;        using the last cetroid to position the mask.  This 
;        can be useful for really noisy data, or data with
;        flat (e.g. saturated) peaks.  Use with caution: it
;        may 'climb' hills and give you multiple hits.
;    quiet: Supress informational messages.
;
; OUTPUTS:
;    f[0,*]: x centroid positions, in pixels.
;    f[1,*[: y centroid positions, in pixels. 
;    f[2,*]: integrated brightness of the features.
;    f[3,*]: square of the radius of gyration
;
; SIDE EFFECTS:
;    Optionally prints diagnostic messages to the terminal.
;
; RESTRICTIONS:
;    This program finds the centers of brightness of bright, circularly
;    symmetric, non-overlapping features on a dark background.  To
;    identify dark features on a bright background, first invert the
;    image.  
;
;    Performance is degraded by non-zero background values,
;    particularly if the background varies across the image.  
;    High-frequency noise also degrades performance.  Both the
;    background and noise can be suppressed by
;    filtering the image with BPASS before running FEATHRE.
;
; PROCEDURE:
;    First, identify the positions of all the local maxima in the
;    image.  Local maxima are the brightest pixels in a circular
;    neighborhood whose diameter is set by DIAMETER.  Within each
;    circular region, the feature's centroid is calculated as
;    the brightness-weighted center of brightness.  The integrated
;    brightness ("mass") and the brightness-weighted radius
;    ("radius of gyration") also are calculated within each circle.
;    If the centroid is found to be more than 0.5 pixels from the
;    original local maximum, the mask can be moved and the centroid
;    recalculated.  This is controlled by the ITERATE keyword.
;    Finally, centroids closer than SEPARATION are merged. 
;
;    This procedure can yield centroid positions with errors of 0.1
;    pixel or better in each dimension for features larger than about
;    5 pixels across.  Achieving this accuracy requires meeting the
;    conditions described in RESTRICTIONS.  Choosing an inappropriate
;    value for DIAMETER or improper settings in BPASS can degrade
;    performance to single-pixel accuracy.
;
;    Finally, FEATURE can select features based on their peak and
;    integrated brightness.  The former is set with the MIN keyword,
;    and the latter with MASSCUT.  Setting these values appropriately
;    greatly increases FEATURE's speed and yields much better
;    rejection of spurious features.  The PICKN keyword also is useful
;    when the number of objects is known.
;
;    Setting DIAMETER, SEPARATION, MIN and MASSCUT is greatly
;    facilitated with the companion program FEATURETOOL.
;
;    The algorithms used in FEATURE and BPASS and a quantitative
;    analysis of their performance are described in
;
;    J.C. Crocker and D.G. Grier, J. Colloid Interface Sci. 179, 298 (1996).
;
; MODIFICATION HISTORY:
;	Original version: feature_stats2: David G. Grier, U of Chicago 1992.
;	Rewritten by John C. Crocker, U of Chicago, optimizing 
;		runtime and measurement error JCC 10/93.
;	Added field keyword JCC 4/94. 
;	Added eccentricity parameter JCC 5/95.
;	Added quiet keyword JCC 12/95.
;	Added iteration, fixed up the radius/diameter fiasco and
;	debugging to improve non-centroid data JCC 4/96.
;       Memory and run-time optimizations DGG 8/99.
;       Added PICKN DGG and E.R. Dufresne 8/99.
;       Fixed error for spheres very near edges DGG 3/2000.
; 11/12/2005. David G. Grier, New York University
;     Major overhaul ...
;     FEATURE now selects each subarray only once -- major speed-up
;     Removed LMX to eliminate redundant code -- another speed-up.
;     Eliminate FOR loop in RSQD.
;     Removed eccentricity, (and thus THETARR) -- rarely used.
;     Updated syntax for IDL 6.X
; 05/20/2006: DGG, improved fidelity to SEPARATION parameter.
;     Overhauled documentation.
; 06/13/2006: DGG: Better respect for QUIET.
;     Don't allow SEP < EXTENT.
; 06/09/2010: DGG. Set COUNT = 0 if no features are found.
;     Fixed bugs with parameter sanity checks.  Documentation fixes.
; 06/10/2010: DGG.  Added COMPILE_OPT statements
; 03/27/2013 DGG More efficient array manipulations.  Revamped message
;   code.  Added usage message.  Remove calls to FIELDOF in favor of
;   array indexing.
;
; Copyright (c) 2006-2013 John C. Crocker, Eric R. Dufresne,
;                           and David G. Grier.
;-

;;;
;
;	RSQD: produce a parabolic mask
;
function rsqd, w, h

COMPILE_OPT IDL2, HIDDEN

if n_params() eq 1 then h = w
xc = float(w-1) / 2.
yc = float(h-1) / 2.
x2 = (findgen(w) - xc)^2
y2 = (findgen(1, h) - yc)^2

return,  rebin(x2, w, h, /sample) + rebin(y2, w, h, /sample)
end

;;;
;
;	FRACSHIFT: barrel shifts a floating point array by a 
;               fractional pixel amount using bilinear interpolation.
;
function fracshift, im, shiftx, shifty

COMPILE_OPT IDL2, HIDDEN

if n_elements(im) eq 1 then return, im

ipx = floor(shiftx)             ; integer part of x shift
ipy = floor(shifty)             ; integer part of y shift
fpx = shiftx - ipx              ; fractional part of x shift
fpy = shifty - ipy              ; fractional part of y shift
if fpx lt 0 then begin
   fpx++ & ipx--                ; make sure fractional parts are positive
endif		
if fpy lt 0 then begin
   fpy++ & ipy--
endif

image = im                      ; preserve input data: use local copy

imagex  = shift(image, ipx+1, ipy  )
imagey  = shift(image, ipx  , ipy+1)
imagexy = shift(image, ipx+1, ipy+1)
image   = shift(image, ipx  , ipy  )

res = ((1. - fpx) * (1. - fpy) * image  ) + $
      ((     fpx) * (1. - fpy) * imagex ) + $
      ((1. - fpx) * (     fpy) * imagey ) + $
      ((     fpx) * (     fpy) * imagexy) 
	
return,res	
end

;;;
;
;	FEATURE:
;
function feature, a, extent, sep,    $
                  min     = min,     $
                  masscut = masscut, $
                  pickn   = pickn,   $
                  field   = field,   $
                  quiet   = quiet,   $
                  iterate = iterate, $
                  maxits  = maxits,  $
                  lmax    = lmax,    $
                  count   = count

COMPILE_OPT IDL2

umsg = 'USAGE: p = feature(image, extent, separation)'

; set flags
quiet = keyword_set(quiet)
field = keyword_set(field)      ; work on a field rather than a frame
iterate = keyword_set(iterate)  ; iterate to improve centroid estimates

; check and process inputs
sz = size(a)
nx = sz[1]                      ; width of image
ny = sz[2]                      ; height of image

extent = floor(extent)          ; diameter of a particle
if (extent mod 2) eq 0 then begin
   extent++
   message, umsg, /inf, noprint = quiet
   message, 'EXTENT must be odd.  Adding 1...', /inf, noprint = quiet
endif

if (n_params() lt 3) then $
  sep = extent + 1              ; separation between particles

if sep le extent then begin
  sep = extent + 1
  message, umsg, /inf, noprint = quiet
  message, 'SEPARATION must be greater than EXTENT: Fixing ...', /inf, noprint = quiet
endif

; derived parameters
radius = float(extent)/2.       ; estimated particle radius
range = floor(sep/2)            ; range over which to search for neighboring features
yrange = (field) ? range/2. : range
yscale = (field) ? 2. : 1.

; numerical options
if n_elements(maxits) ne 1 then maxits = 10 ; maximum number of iterations
if n_elements(masscut) ne 1 then masscut = 0 ; minimum acceptable integrated brightness

if ~keyword_set(min) then begin ; minimum acceptable pixel intensity
   h = histogram(a)
   goal = 0.64 * total(h)
   min = 0
   val = h[min]
   while val le goal do begin
      min++
      val += h[min]
   endwhile
   message, "Setting MIN to " + strcompress(min), /inf, noprint = quiet
endif

; the array of results
res = fltarr(4)

; find local maxima
mmask = rsqd(sep) lt range^2
if field then mmask = mmask[*, 1:*:2] ; odd field by default
b = byte(a)
c = dilate(b, mmask, /gray)
r = where((b eq c) and (b ge min), count)
if count lt 1 then begin
   message, "No local maxima were brighter than MIN", /inf, noprint = quiet
   return, res
endif

; local maxima provide initial estimates for particle positions
x  = float(r mod nx)
y  = float(floor(r / nx)) 

; some local maxima will be too close to the edge -- eliminate them
good = where((x ge range) and (x lt (nx-range)) and $
             (y ge yrange) and (y lt (ny-yrange)), lmax)
if lmax lt 1 then begin
   message, "All local maxima were too close to edge", /inf, noprint = quiet
   count = 0
   return, res
endif
x = x[good]
y = y[good]
message, strcompress(lmax) + ' local maxima found.', /inf, noprint = quiet

; corners of regions around each local maximum
xl = x - floor(extent/2) 
xh = xl + extent - 1
yl = (field) ? y - floor(extent/4) : y - floor(extent/2) 
yh = (field) ? yl + floor(extent/2) - 1 : yl + extent - 1

; set up some masks
rsq = rsqd(extent)
mask = rsq le radius^2
xmask = rebin(findgen(extent) - (extent - 1)/2., extent, extent, /sample)
ymask = transpose(xmask)
if field then begin
   mask = mask[*, 1:*:2]
   xmask = xmask[*, 1:*:2]
   ymask = ymask[*, 1:*:2]
   rsq = rsq[*, 1:*:2]
endif
xmask *= mask
ymask *= mask
rmask = rsq * mask + 1./6.

res = fltarr(4, lmax, /nozero)

; calculate sub-pixel centroid and other characteristics for each feature
for i = 0L, lmax-1L do begin
   xi =  x[i] & yi = y[i]
   xli = xl[i] & xhi = xh[i]
   yli = yl[i] & yhi = yh[i]

   suba = a[xli:xhi, yli:yhi]   ; sub-image around candidate feature
   m = total(suba * mask)       ; integrated brightness

   if m lt masscut then continue ; too small: m = 0 and rg = 0 at loop's end
                                ; doing the test here eliminates
                                ; calculations for spurious features.

                                ; displacement of centroid from center of mask
    xc = total(suba * xmask) / m
    yc = total(suba * ymask) / (m * yscale)
        
    if iterate then begin       ; iterate to improve position estimate
       its = 0
       repeat begin
          shifted = 0
          if abs(xc) gt 0.6 then begin
             shifted = 1
             dx = round(xc)
             xli = (xli + dx) > 0
             xhi = (xhi + dx) <  (nx-1) > (xli+1)
             xi = (xi + dx) > xli <  xhi
          endif
          if abs(yc) gt 0.6 then begin
             shifted = 1
             dy = round(yc)
             yli = (yli + dy) > 0
             yhi = (yhi + dy) <  (ny-1) > (yli+1)
             yi = (yi + dy) > yli <  yhi
          endif
          if shifted then begin     
             suba = a[xli:xhi, yli:yhi]
             m = total(suba * mask)
             xc = total(suba * xmask) / m
             yc = total(suba * ymask) / (m * yscale)
            endif
          its++
       endrep until not shifted or (its eq maxits)
    endif

    rg = total(suba * rmask) / m ; radius of gyration
    res[*, i] = [xi+xc, (yi+yc)*yscale, m, rg]
endfor

; keep only non-trivial features
w = where(res[2, *] gt 0, count)
if count gt 0 then $
   res = res[*, w] $
else $
   return, res

; some features might have converged to the same point
; eliminate duplicates
hash =  floor(res[0, *]/sep) + nx * floor(res[1, *]/sep)
ndx = uniq(hash, sort(hash))
count = n_elements(ndx)
res = res[*, ndx]
message, strcompress(count) + " unique features above threshold", /inf, noprint = quiet

; select the brightest features
if keyword_set(pickn) then begin
   if count lt pickn then $
      message, "PICKN: Ignored: Fewer than " + strtrim(pickn,2) + $
               " features to choose from.", /inf, noprint = quiet $
    else begin
      order = sort(res[2, *])   ; sort by integrated brightness
      good = order[count-pickn:*]    
      res = res[*, good]
      message, strcompress(pickn) + " features selected", /inf, noprint = quiet
   endelse
endif

return, res
end
