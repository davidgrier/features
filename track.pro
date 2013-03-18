; ********* start of track.pro
;
; see http://www.physics.emory.edu/~weeks/idl
;   for more information
;
;+
; NAME:
;    track
;
; PURPOSE:
;    Constructs n-dimensional trajectories from a scrambled list of
;    particle coordinates determined at discrete times (e.g. in
;    consecutive video frames).
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    result = track(positionlist, maxdisp [, inipos=inipos, 
;                   memory=memory, dim=dim, verbose=verbose, 
;                   goodenough=goodenough, quiet=quiet ])
; INPUTS:
;    positionlist: an array listing the scrambled coordinates and data 
;        of the different particles at different times, such that:
;        positionlist[0:d-1,*]: contains the d coordinates and
;            data for all the particles, at the different times.
;        positionlist[d,*]: contains the time t that the position 
;            was determined (e.g. frame number.  These values must 
;            be monotonically increasing and uniformly gridded in time.
;    maxdisp: an estimate of the maximum distance that a particle 
;        would move in a single time interval.(see Restrictions)
;
; KEYWORD PARAMETERS
;    inipos: if the user wants to track only a subset of the 
;        particles, this argument is set to an array (d,n)
;        which contains the d dimensional initial positions of
;        the n particles to be tracked. Other 'new' particles
;        will then NOT be added.
;    memory: this is the number of time steps that a particle can be
;        'lost' and then recovered again.  If the particle reappears
;        after this number of frames has elapsed, it will be
;        tracked as a new particle. The default setting is zero.
;        this is useful if particles occasionally 'drop out' of
;        the data.
;    dim: if the user would like to unscramble non-coordinate data
;        for the particles (e.g. apparent radius of gyration for
;        the particle images), then positionlist should
;        contain the position data in positionlist[0:dim-1,*]
;        and the extra data in positionlist[dim:d-1,*]. It is then
;        necessary to set dim equal to the dimensionality of the
;        coordinate data to so that the track knows to ignore the
;        non-coordinate data in the construction of the 
;        trajectories. The default value is two.
;    goodenough: set this keyword to eliminate all trajectories with
;        fewer than goodenough valid positions.  This is useful
;        for eliminating very short, mostly 'lost' trajectories
;        due to blinking 'noise' particles in the data stream.
;
; KEYWORD FLAGS
;    verbose: set this keyword for more informational messages.
;
;    quiet:  use /quiet to not print any messages
;
; OUTPUTS:
;    result: a list containing the original data rows sorted 
;        into a series of trajectories.  To the original input 
;        data structure there is appended an additional column 
;        containing a unique 'id number' for each identified 
;        particle trajectory.  The result array is sorted so 
;        rows with corresponding id numbers are in contiguous 
;        blocks, with the time variable a monotonically
;        increasing function inside each block.  For example:
;		
;		For the input data structure (positionlist):
;			(x)	     (y)	  (t)
;		pos = 3.60000      5.00000      0.00000
;		      15.1000      22.6000      0.00000
;		      4.10000      5.50000      1.00000	
;		      15.9000      20.7000      2.00000
;		      6.20000      4.30000      2.00000
;
;		IDL> res = track(pos,5,mem=2)
;
;		track will return the result 'res'
;			(x)	     (y)	  (t) 	       (id)
;		res = 3.60000      5.00000      0.00000      0.00000
;		      4.10000      5.50000      1.00000      0.00000
;		      6.20000      4.30000      2.00000      0.00000
;		      15.1000      22.6000      0.00000      1.00000
;		      15.9000      20.7000      2.00000      1.00000
;
;        NB: for t=1 in the example above, one particle temporarily
;        vanished.  As a result, the trajectory id=1 has one time
;        missing, i.e. particle loss can cause time gaps to occur 
;        in the corresponding trajectory list. In contrast:
;
;		IDL> res = track(pos,5)
;
;		track will return the result 'res'
;			(x)	     (y)	  (t) 	       (id)
;		res = 15.1000      22.6000      0.00000      0.00000
;      		      3.60000      5.00000      0.00000      1.00000
;  		      4.10000      5.50000      1.00000      1.00000
; 		      6.20000      4.30000      2.00000      1.00000
; 		      15.9000      20.7000      2.00000      2.00000
;	
;        where the reappeared 'particle' will be labelled as new
;        rather than as a continuation of an old particle since
;        mem=0.  It is up to the user to decide what setting of 
;        'mem' will yeild the highest fidelity tracking.
;	
; SIDE EFFECTS:
;    Produces informational messages.  Can be memory intensive for
;    extremely large data sets.
; RESTRICTIONS:
;    maxdisp should be set to a value somewhat less than the mean 
;    spacing between the particles. As maxdisp approaches the mean
;    spacing the runtime will increase significantly. The function 
;    will produce an error message: "Excessive Combinatorics!" if
;    the run time would be too long, and the user should respond 
;    by re-executing the function with a smaller value of maxdisp.
;    Obviously, if the particles being tracked are frequently moving
;    as much as their mean separation in a single time step, this
;    function will not return acceptable trajectories.
; PROCEDURE:
;    Given the positions for n particles at time t(i), and m possible
;    new positions at time t(i+1), this function considers all possible 
;    identifications of the n old positions with the m new positions,
;    and chooses that identification which results in the minimal total
;    squared displacement. Those identifications which don't associate
;    a new position within maxdisp of an old position ( particle loss )
;    penalize the total squared displacement by maxdisp^2. For non-
;    interacting Brownian particles with the same diffusivity, this
;    algorithm will produce the most probable set of identifications 
;    ( provided maxdisp >> RMS displacement between frames ).
;    In practice it works reasonably well for systems with oscillatory,
;    ballistic, correlated and random hopping motion, so long as single 
;    time step displacements are reasonably small.  NB: multidimensional
;    functionality is intended to facilitate tracking when additional
;    information regarding target identity is available (e.g. size or 
;    color).  At present, this information should be rescaled by the
;    user to have a comparable or smaller (measurement) variance than 
;    the spatial displacements.
;
; MODIFICATION HISTORY:
; 02/93 Written by John C. Crocker, University of Chicago (JFI).
; 07/93 JCC fixed bug causing particle loss and improved performance
;    for large numbers of (>100) particles.
; 11/93 JCC improved speed and memory performance for large
;    numbers of (>1000) particles (added subnetwork code).
; 03/94 JCC optimized run time for trivial bonds and d<7. (Added
;    d-dimensional raster metric code.)
; 08/94 JCC added functionality to unscramble non-position data
;    along with position data.
; 09/94 JCC rewrote subnetwork code and wrote new, more efficient 
;    permutation code.
; 05/95 JCC debugged subnetwork and excessive combinatorics code.
; 12/95 JCC added memory keyword, and enabled the tracking of
;    newly appeared particles.
; 03/96 JCC made inipos a keyword, and disabled the adding of 'new'
;    particles when inipos was set.
; 03/97 JCC added 'add' keyword, since Chicago users didn't like 
;    having particle addition be the default. 
; 09/97 JCC added 'goodenough' keyword to improve memory efficiency
;    when using the 'add' keyword and to filter out bad tracks.
; 10/97 JCC streamlined data structure to speed runtime for >200 
;    timesteps.  Changed 'quiet' keyword to 'verbose'. Made
;    time labelling more flexible (uniform and sorted is ok).
; 09/98 JCC switched trajectory data structure to a 'list' form,
;    resolving memory issue for large, noisy datasets.
; 09/17/1998 Eric Weeks, Emory University, luberize code.
; 02/99 JCC added Eric Weeks's 'uberize' code to post-facto 
;    rationalize the particle id numbers, removed 'add' keyword.
;    First public release of track.pro
; 03/30/2010 David G. Grier, New York University: Modernized array
;    notation.  Small code modernizations.  Formatting.
;    Moved luberize code into main procedure.  Replaced UNQ with
;    IDL system routine.
; 05/25/2012 Eric Weeks: Updated documentation.  Change findgens to
;    lindgens to accommodate very large data sets.  Long counters in
;    for loops.  Added QUIET keyword.
; 03/17/2013 DGG Added COMPILE_OPT.  Simplify analysis of time steps.
;    Replace bitwise AND, OR and NOT operators with &&, || and ~.
;    Replace keyword_set() with numeric tests, where appropriate.
;    Simplify "goodenough" tests.  Simplify "inipos" tests.
;    Calculate CUBE in "the IDL way".
;
;	This code 'track.pro' is copyright 1999, by John C. Crocker. 
;	It should be considered 'freeware'- and may be distributed freely 
;	(outside of the military-industrial complex) in its original form 
;	when properly attributed.
;
;-
function track, xyzs, maxdisp, $
                inipos = inipos, $
                memory = memory, $
                dim = dim,$
                goodenough = goodenough, $
                verbose = verbose, $
                quiet = quiet

COMPILE_OPT IDL2

dd = n_elements(xyzs[*,0]) - 1
if ~isa(dim, /number, /scalar) then begin
   dim = 2 < dd
   message,' Setting dim = ' + strtrim(dim,2) + ' by default', /inf
endif

if ~isa(memory, /number, /scalar) then memory = 0

; check the input time vector is ok, i.e. sorted and uniform
t = reform(xyzs[dd,*])
dt = t[1:*] - t
if total(dt lt 0) ne 0 then $
   message, ' Error- Times are not monotonically increasing!'

w = where(dt gt 0, nsteps)
if nsteps eq 0 then $
   message, ' Error- All data have the same time stamp!'

if ~array_equal(dt[w], dt[w[0]]) then $
   message, ' Warning- Times are not evenly gridded!', /inf

nsteps++                        ; number of time steps

doprune = isa(goodenough, /number, /scalar) && (goodenough gt 0)
doinipos = isa(inipos, /number, /array)

; partition the input data by unique times
res = uniq(t)
res = [0, res+1, n_elements(t)]

; get the initial positions
ngood = res[1] - res[0]
eyes = lindgen(ngood) + res[0]

if doinipos then begin
   pos = inipos[0:dim-1, *]
   istart = 0L 
   n = n_elements(pos[0, *]) 
endif else begin
   pos = xyzs[0:dim-1, eyes]
   istart = 1L                  ; we don't need to track t=0.
   n = ngood
endelse

;how long are the 'working' copies of the data?
zspan = 50
if n gt 200 then zspan = 20
if n gt 500 then zspan = 10

resx = lonarr(n, zspan) - 1
bigresx = lonarr(n, nsteps) - 1
mem = lonarr(n)
uniqid = findgen(n)
maxid = n
olist = [0., 0.]

if doprune then begin
   dumphash = bytarr(n)
   nvalid = intarr(n)
endif

; we may not need to track the first time step!
if ~doinipos then begin
   resx[*, 0] = eyes
   if doprune then nvalid++
endif

; set up some nice constants
maxdisq = maxdisp^2
verbose = keyword_set(verbose)
quiet = keyword_set(quiet)

;Use fancy code for large n, small d
notnsqrd = (sqrt(n*ngood) ge 200) && (dim lt 7)

if notnsqrd then begin
;   construct the vertices of a 3x3x3... d-dimensional hypercube
   cube = rebin(lindgen(1, 3^dim), dim, 3^dim)
   for j = 1, dim-1 do cube[j, *] /= 3^j
   cube mod= 3

;   calculate a blocksize which may be greater than maxdisp, but which
;   keeps nblocks reasonably small.   
   volume = 1
   for d = 0, dim-1 do begin
      minn = min(xyzs[d, w], max = maxx)
      volume *= maxx - minn
   endfor
   blocksize = maxdisp > (float(volume)/(20.*ngood))^(1./dim)
endif

;   Start the main loop over the frames.
for i = istart, nsteps-1 do begin

   ispan = i mod zspan   
      
;   Get the new particle positions.
   m = res[i+1] - res[i]
   eyes = lindgen(m) + res[i]

   if m gt 0 then begin
      xyi = xyzs[0:dim-1, eyes]
      found = bytarr(m)

;   THE TRIVIAL BOND CODE BEGINS   
      if notnsqrd then begin   
;   Use the raster metric code to do trivial bonds

;   construct "s", a one dimensional parameterization of the space 
;   ( which consists of the d-dimensional raster scan of the volume.)
         abi = long(xyi/blocksize)
         abpos = long(pos/blocksize)
         si = lonarr(m)
         spos = lonarr(n)
         dimm = lonarr(dim, /nozero)
         coff = 1L
         for j = 0, dim-1 do begin
            minn = min([[abi[j,*]],[abpos[j,*]]], max = maxx)
            abi[j,*] -= minn
            abpos[j,*] -= minn
            dimm[j] = maxx - minn + 1
            si += abi[j,*]*coff
            spos += abpos[j,*]*coff
            coff *= dimm[j]
         endfor
         nblocks = coff         ; the # of blocks in the volume
      
;   trim down (intersect) the hypercube if its too big to fit in the
;   particle volume. (i.e. if dimm(j) lt 3)
         cub = cube
         deg = where(dimm lt 3, ninside)
         if ninside gt 0 then begin
            for j = 0, ninside-1 do $
               cub = cub[*, where(cub[deg[j],*] lt dimm[deg[j]])]
         endif
      
;   calculate the "s" coordinates of hypercube (with a corner @ the origin)
         scube = lonarr(n_elements(cub[0,*]))
         coff = 1L
         for j = 0, dim-1 do begin
            scube += cub[j,*]*coff
            coff *= dimm[j]
         endfor   
      
;   shift the hypercube "s" coordinates to be centered around the origin
         coff = 1L
         for j = 0, dim-1 do begin
            if dimm[j] ge 3 then $
               scube -= coff
            coff *= dimm[j]
         endfor   
         scube = (scube + nblocks) mod nblocks         
      
;   get the sorting for the particles by their "s" positions.
         isort = sort(si)
      
;   make a hash table which will allow us to know which new particles
;   are at a given si.
         strt = intarr(nblocks) - 1
         fnsh = intarr(nblocks, /nozero)
         for j = 0, m-1 do begin
            ndx = si[isort[j]]
            if strt[ndx] eq -1 then $
               strt[ndx] = j
            fnsh[ndx] = j
         endfor      

;   loop over the old particles, and find those new particles in the 'cube'.
         coltot = intarr(m)
         rowtot = intarr(n)
         which1 = intarr(n, /nozero)
         for j = 0L, n-1L do begin
            map = fix(-1)
            s = (scube + spos[j]) mod nblocks
            w = where(strt[s] ne -1, ngood)
            if ngood ne 0 then begin
               s = s[w]
               for k = 0L, ngood-1L do $
                  map = [map,isort[strt[s[k]]:fnsh[s[k]]]] 
               map = map[1:*]
               
               ; find those trivial bonds
               distq = fltarr(n_elements(map))
               for d = 0, dim-1 do $
                  distq += (xyi[d,map] - pos[d,j])^2  
               ltmax = distq lt maxdisq
               rowtot[j] = total(ltmax)               

               if rowtot[j] ge 1 then begin
                  w = where(ltmax)
                  coltot[map[w]] += 1 
                  which1[j] = map[w[0]]
               endif   
            endif   
         endfor

         ntrk = fix(n - total(rowtot eq 0))
         w = where(rowtot eq 1, ngood)
         if ngood ne 0 then begin
            ww = where(coltot[which1[w]] eq 1, ngood)
            if ngood ne 0 then begin
               ndx = w[ww]
               resx[ndx,ispan] = eyes[which1[ndx]]
               found[which1[ndx]] = 1B
               rowtot[ndx] = 0
               coltot[which1[ndx]] = 0
            endif
         endif
         labely = where(rowtot gt 0, ngood)
         if ngood ne 0 then begin
            labelx = where(coltot gt 0)
            nontrivial = 1
         endif else $
            nontrivial = 0
         
      endif else begin
;   or: Use simple N^2 time routine to calculate trivial bonds      

         ; let's try a nice, loopless way!
         ; don't bother tracking permanently lost guys.
         wh = where(pos[0,*] ge 0, ntrack)
         if ntrack eq 0 then $
            message, 'Warning: No valid particles to track!'
         
         xmat = make_array(m, ntrack, /long, /index) mod m
         ymat = transpose(make_array(ntrack, m, /long, /index) mod ntrack)
         for d = 0, dim-1 do begin
            x = reform(xyi[d,*])
            y = reform(pos[d,wh])
            if d eq 0 then $
               dq = (x[xmat] - y[ymat])^2 $
            else $
               dq += (x[xmat] - y[ymat])^2
         endfor

         ltmax = dq lt maxdisq

         ; figure out which trivial bonds go with which
         rowtot = intarr(n)
         rowtot[wh] = total(ltmax, 1)
         coltot = (ntrack gt 1) ? total(ltmax, 2) : ltmax
         which1 = intarr(n)
         for j = 0L, ntrack-1L do begin
            mx = max(ltmax[*,j], w) ; max is faster than where
            which1[wh[j]] = w
         endfor

         ntrk = fix(n - total(rowtot eq 0))
         w = where(rowtot eq 1, ngood)
         if ngood ne 0 then begin
            ww = where(coltot[which1[w]] eq 1, ngood)
            ndx = w[ww]
            if ngood ne 0 then begin
               resx[ndx,ispan] = eyes[which1[ndx]]
               found[which1[ndx]] = 1B
               rowtot[ndx] = 0
               coltot[which1[ndx]] = 0
            endif
         endif

         labely = where(rowtot gt 0, ngood)
         if ngood ne 0 then begin
            labelx = where(coltot gt 0)
            nontrivial = 1
         endif else $
            nontrivial = 0
         
      endelse   
;   THE TRIVIAL BOND CODE ENDS         
      if nontrivial then begin
      
         xdim = n_elements(labelx)
         ydim = n_elements(labely)
         
;   make a list of the non-trivial bonds            
         bonds = lonarr(2)
         bondlen = 0.   
         for j = 0, ydim-1 do begin
            distq = fltarr( xdim )
            for d = 0, dim-1 do begin
               distq += (xyi[d,labelx] - pos[d,labely[j]])^2  
            endfor
            w = where(distq lt maxdisq, ngood)
            bonds = [[bonds], [transpose(w), lonarr(1,ngood)+j]]
            bondlen = [bondlen, distq[w]]   
         endfor   
         bonds = bonds[*,1:*]
         bondlen = bondlen[1:*]
         numbonds = n_elements(bonds[0,*])
         mbonds = bonds      
         
         if (xdim > ydim) LT 4 then begin
            nclust = 1
            maxsz = 0
            mxsz = xdim
            mysz = ydim
            bmap = lonarr(n_elements(bonds[0,*])) - 1
         endif else begin
;   THE SUBNETWORK CODE BEGINS            

            lista = intarr(numbonds)
            listb = intarr(numbonds)
            nclust = 0
            maxsz = 0
            thru = xdim
      
            while thru NE 0 do begin               
;   the following code extracts connected sub-networks of the non-trivial 
;   bonds.  NB: lista/b can have redundant entries due to 
;   multiple-connected subnetworks.
               w = where(bonds[1,*] GE 0)
               lista[0] = bonds[1,w[0]]
               listb[0] = bonds[0,w[0]]
               bonds[*,w[0]] = -(nclust+1)
               adda  = 1 & addb  = 1
               donea = 0 & doneb = 0
               repeat begin
                  if (donea NE adda) then begin
                     w = where(bonds[1,*] EQ lista[donea], ngood)   
                     if ngood NE 0 then begin
                        listb[addb:addb+ngood-1] = bonds[0,w]
                        bonds[*,w] = -(nclust+1)
                        addb += ngood
                     endif
                     donea += 1
                  endif
                  if (doneb NE addb) then begin
                     w = where(bonds[0,*] EQ listb[doneb], ngood)   
                     if ngood NE 0 then begin
                        lista[adda:adda+ngood-1] = bonds[1,w]
                        bonds[*,w] = -(nclust+1)
                        adda += ngood
                     endif
                     doneb += 1
                  endif
               endrep until (donea EQ adda) && (doneb EQ addb)
               ; a thing of beauty is a joy forever.

               xsz = n_elements(uniq(listb[0:doneb-1],sort(listb[0:doneb-1])))
               ysz = n_elements(uniq(lista[0:donea-1],sort(lista[0:donea-1])))
               
               if xsz*ysz GT maxsz then begin
                  maxsz = xsz*ysz
                  mxsz = xsz
                  mysz = ysz
               endif
               
               thru -= xsz
               nclust += 1
            
            endwhile
            bmap = reform(bonds[0,*])                           
         endelse            
;   THE SUBNETWORK CODE ENDS
         if verbose then begin
            message, strcompress(i)+': '+'Permuting'+ $
                     strcompress(nclust)+' network'+ $
                     (nclust GT 1) ? 's' : ' ', /inf 
            message,'      Max. network:'+strcompress(mxsz)+' x'+$
                    strcompress(mysz), /inf
            if keyword_set(add) then $
               message,'      Tracking'+ $
                       strcompress(n)+' particles.', /inf
         endif            
            
;   THE PERMUTATION CODE BEGINS
         for nc = 0, nclust-1 do begin

            w = where(bmap EQ -(nc+1), nbonds)
            bonds = mbonds[*, w]
            lensq = bondlen[w]
               
            uold = bonds[0,uniq(bonds[0,*],sort(bonds[0,*]))]
            nold = n_elements(uold)
            unew = bonds[1,uniq(bonds[1,*])]
            nnew = n_elements(unew)
            
            ; check that runtime is not excessive
            if nnew gt 5 then begin
               rnsteps = 1D
               for ii = 0, nnew-1 do begin
                  rnsteps *= n_elements(where(bonds[1,*] eq unew[ii]))
                  if rnsteps gt 5.e4 then $
                     message, $
                     ' Warning: difficult combinatorics encountered.', /inf
                  if rnsteps gt 2.e5 then $
                     message, $
                     ' Excessive Combinatorics! Try reducing maxdisp.'  
               endfor   
            endif
            
            st = intarr(nnew) & fi = intarr(nnew)            
            h = intarr(nbonds) & ok = intarr(nold) + 1
            nlost = (nnew - nold) > 0               
               
            for ii = 0, nold-1 do $
               h[where(bonds[0,*] EQ uold[ii])] = ii
            
            st[0] = 0 & fi[nnew-1] = nbonds-1
            if nnew gt 1 then begin
               sb = reform(bonds[1,*])
               sbr = shift(sb, 1)
               sbl = shift(sb, -1)
               st[1:*] = where(sb[1:*] ne sbr[1:*]) + 1
               fi[0:nnew-2] = where(sb[0:nbonds-2] ne sbl[0:nbonds-2])
            endif
               
            checkflag = 0
            repeat begin               
               pt = st - 1
               lost = intarr(nnew) 
               who = 0 & losttot = 0
               mndisq = nnew*maxdisq
               
               repeat begin
                  if pt[who] ne fi[who] then begin
                     w = where(ok[h[pt[who]+1:fi[who]]], ngood)
                     if ngood gt 0 then begin
                        if pt[who] ne st[who] - 1 then $
                           ok[h[pt[who]]] = 1
                        pt[who] += 1 + w[0]
                        ok[h[pt[who]]] = 0
                        if who eq nnew-1 then begin
                           ; place #1 calc. tot. sqr. disp.
                           ww = where(lost eq 0)
                           dsq = total(lensq[pt[ww]]) + losttot*maxdisq   

                           if dsq lt mndisq then begin
                              minbonds = pt[ww]
                              mndisq = dsq
                           endif
                        endif else $
                           who += 1
                     endif else begin
                        if ~lost[who] && (losttot ne nlost) then begin
                           lost[who] = 1
                           losttot += 1
                           if pt[who] ne st[who] - 1 then $
                              ok[h[pt[who]]] = 1       
                           if who eq nnew-1 then begin
                           ; place #2 calc. tot. sqr. disp.
                              ww = where(lost eq 0)
                              dsq = total(lensq[pt[ww]]) + losttot*maxdisq   
                              
                              if dsq lt mndisq then begin
                                 minbonds = pt[ww]
                                 mndisq = dsq
                              endif
                           endif else $
                              who += 1
                           ; Fight the 'c' power- long live IDL!
                        endif else begin
                           if pt[who] ne st[who] - 1 then $
                              ok[h[pt[who]]] = 1
                           pt[who] = st[who] - 1
                           if lost[who] then begin
                              lost[who] = 0
                              losttot -= 1
                           endif
                           who -= 1
                        endelse
                     endelse   
                  endif else begin
                     if ~lost[who] && (losttot ne nlost) then begin
                        lost[who] = 1
                        losttot += 1
                        if pt[who] ne st[who] - 1 then $
                           ok[h[pt[who]]] = 1   
                        if who eq nnew-1 then begin
                           ; place #3 calc. tot. sqr. disp.
                           ww = where(lost eq 0)
                           dsq = total(lensq[pt[ww]]) + losttot*maxdisq   
                           
                           if dsq lt mndisq then begin
                              minbonds = pt[ww]
                              mndisq = dsq
                           endif
                        endif else $
                           who += 1
                     endif else begin
                        if pt[who] ne st[who] - 1 then $
                           ok[h[pt[who]]] = 1
                        pt[who] = st[who] -1
                        if lost[who] then begin
                           lost[who] = 0
                           losttot -= 1
                        endif
                        who -= 1
                     endelse
                  endelse         
               endrep until who eq -1 
               
               checkflag += 1
               if checkflag eq 1 then begin
;   we need to check that our constraint on nlost is not forcing us away from 
;   the minimum id's
                  plost = fix(mndisq/maxdisq) < nnew-1
                  if plost gt nlost then $
                     nlost = plost $
                  else $
                     checkflag = 2
               endif
               
            endrep until checkflag eq 2
            
;   update resx using the minimum bond configuration               
            resx[labely[bonds[1,minbonds]],ispan] = $
               eyes[labelx[bonds[0,minbonds]]]
            found[labelx[bonds[0,minbonds]]] = 1
         endfor
;   THE PERMUTATION CODE ENDS
         
      endif else if verbose then $
         message, strcompress(i)+': Only trivial networks', /inf

;     here we want to update our initial position estimates
      w = where(resx[*,ispan] ge 0, nww)
      if nww gt 0 then begin
         pos[*,w] = xyzs[0:dim-1,resx[w,ispan]]
         if doprune then nvalid[w]++
      endif else $
         message, ' Warning, tracking zero particles!', /inf

;     we need to add new guys, as appropriate.
      newguys = where(found eq 0, nnew)
      if (nnew gt 0) && ~doinipos then begin
         newarr = fltarr(nnew, zspan) - 1.
         resx = [resx,newarr]
         resx[n:*,ispan] = eyes[newguys]
         pos = [[pos],[xyzs[0:dim-1,eyes[newguys]]]]
         mem = [mem,bytarr(nnew)]
	 uniqid = [uniqid,findgen(nnew)+maxid]
	 maxid += nnew
         if doprune then begin
            dumphash = [dumphash,bytarr(nnew)]
            nvalid = [nvalid,intarr(nnew)+1]
         endif
         n += nnew
      endif

   endif else $
      message, ' Warning- No positions found for t='+strcompress(i)+"!", /inf
   
;   update the 'memory' array
   w = where(resx[*,ispan] ne -1, nok)
   if nok ne 0 then $
      mem[w] = 0                ; guys get reset if they're found
   mem += resx[*,ispan] eq -1

;  if a guy has been lost for more than memory times, mark him as permanently
;  lost.  For now, set these guys to pos = ( -maxdisp, -maxdisp, ... ),
;  so we can never track them again. It would be better to make a smaller
;  pos, but then we'd have to change 'n', which would be gnarly.
   wlost = where(mem eq memory+1, nlost)
   if nlost gt 0 then begin
      pos[*,wlost] = -maxdisp
      ; check to see if we should 'dump' newly lost guys
      if doprune then begin
         wdump = where(nvalid[wlost] lt goodenough, ndump)
         if ndump gt 0 then $
            dumphash[wlost[wdump]] = 1B
      endif
   endif

;  we need to insert the working copy of resx into the big copy bigresx
;  do our house keeping every zspan time steps (dumping bad lost guys)
   if (ispan eq zspan-1) or (i eq nsteps-1) then begin

;  if a permanently lost guy has fewer than goodenough valid positions
;  then we 'dump' him out of the data structure- this largely alleviates
;  memory problems associated with the 'add' keyword and 'noise' particles
;  To improve speed- do it infrequently.
      ; in case we've added some we need to pad out bigresx too
      nold = n_elements(bigresx[*,0])
      nnew = n - nold
      if nnew gt 0 then begin   
         newarr = fltarr(nnew, nsteps) - 1.
         bigresx = [bigresx,newarr]
      endif
      if doprune then begin
         if (total(dumphash) gt 0) then begin
            if verbose then $
               message, 'Dumping bad trajectories...', /inf
            wkeep = where(dumphash eq 0, nkeep)
            resx = resx[wkeep,*]
            bigresx = bigresx[wkeep,*] ; this really hurts runtime
            pos = pos[*,wkeep]
            mem = mem[wkeep]
            uniqid = uniqid[wkeep]
            nvalid = nvalid[wkeep]
            n = nkeep
            dumphash = bytarr(nkeep)
         endif
      endif
      if ~verbose && ~quiet then $
         message, strcompress(i+1)+' of'+strcompress(nsteps)+' done. '+ $
                  'Tracking'+strcompress(ntrk)+' particles, '+ $
                  strcompress(n)+' tracks total.', /inf
      bigresx[*,i-ispan:i] = resx[*,0:ispan]
      resx = lonarr(n,zspan) - 1

;  We should pull permanently lost guys, parse them and concat them
;  onto the 'output list', along with their 'unique id' number to
;  make scanning the data files a little easier.  Do infrequently.
      wpull = where(pos[0,*] eq -maxdisp, npull)
      if npull gt 0 then begin
      ;   pos( 0, wpull ) = -2*maxdisp
         lillist = [0.,0.]
         for ipull = 0, npull-1 do begin
            wpull2 = where(bigresx[wpull[ipull],*] ne -1, npull2)
            lillist = [[lillist], $
                       [bigresx[wpull[ipull],wpull2], $
                        fltarr(1,npull2)+uniqid[wpull[ipull]]]]
         endfor
         olist = [[olist],[lillist[*,1:*]]]
      endif
;     now get rid of the guys we don't need any more....  
;     but watch out for when we have no valid particles to track!    
      wkeep = where(pos[0,*] ge 0, nkeep)
      if nkeep eq 0 then $
         print," We're going to crash now, no particles...."
      resx = resx[wkeep,*]
      bigresx = bigresx[wkeep,*] ; this really hurts runtime
      pos = pos[*,wkeep]
      mem = mem[wkeep]
      uniqid = uniqid[wkeep]
      n = nkeep
      dumphash = bytarr(nkeep)
      if doprune then $
         nvalid = nvalid[wkeep]
   endif

endfor  ; the big loop over nsteps time steps....

;  make a final scan for short trajectories that weren't lost at the end. 
if doprune then begin
   nvalid = total(bigresx ge 0, 2)
   wkeep = where(nvalid ge goodenough, nkeep)
   if nkeep lt n then begin
       bigresx = bigresx[wkeep,*]
       n = nkeep
       uniqid = uniqid[wkeep]
       pos = pos[*,wkeep]
   endif
endif

;  make the final scan to 'pull' everybody else into the olist.
wpull = where(pos[0,*] ne -2*maxdisp, npull)
if npull gt 0 then begin
   lillist = [0.,0.]
   for ipull = 0, npull-1 do begin
      wpull2 = where(bigresx[wpull[ipull],*] ne -1, npull2)
      lillist = [[lillist], $
                 [bigresx[wpull[ipull],wpull2], $             
                  fltarr(1,npull2)+uniqid[wpull[ipull]]]]
   endfor
   olist = [[olist],[lillist[*,1:*]]]
endif
olist = olist[*,1:*]

;  free up a little memory for the final step!
bigresx = 0
resx = 0

; need to make up a result array!
if verbose then $
   message, 'Preparing result array...', /inf
nolist = n_elements(olist[0,*])
res = fltarr(dd+2, nolist)
for j = 0, dd do begin
   res[j,*] = xyzs[j,olist[0,*]]
endfor
res[dd+1,*] = olist[1,*]

ndat = n_elements(res[*,0]) - 1
u = uniq(res[ndat,*])
ntracks = n_elements(u)
u = [-1,u]
for i = 1L, ntracks do $
   res[ndat,u[i-1]+1:u[i]] = i-1

return, res
end
; ***************** end of track.pro

