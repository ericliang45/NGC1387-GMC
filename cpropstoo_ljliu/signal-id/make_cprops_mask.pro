pro make_cprops_mask $
   , infile = infile $
   , indata = cube $
   , inmask = inmask $
   , outmask = mask $
   , outfile=outfile $   
   , outisland=outisland $   
   , outmomzero=outmomzero $   
   , outmomone=outmomone $   
   , rmsfile = rmsfile $
   , inrms = rms $
;  CONDITIONS FOR THE MASK
   , prior = prior $
   , hi_thresh = hi_thresh $
   , hi_nchan = hi_nchan $
   , min_pix = min_pix $
   , min_area = min_area $   
;  EXTEND TO INCLUDE FAINTER ENVELOPES
   , lo_thresh = lo_thresh $
   , lo_nchan = lo_nchan $
;  DILATION
   , grow_xy = grow_xy $
   , grow_z = grow_z $
;  CLIP REGIONS BASED ON RMS LEVEL
   , clip_rms = clip_rms $
;  INVERT NEGATIVE->POSITIVE (USEFUL TO CHECK FALSE POSITIVES)
   , invert = invert $
;  GIVE OUTPUT
   , verbose = verbose $
;  CLIP EDGE
   , edge_pc = edge_pc $
   , min_edge_area = min_edge_area $
   , indist = indist $
;  CONVERT MASK TO FLOAT BEFORE OUTPUT (CASA WONT READ BYTE ARRAYS)
   , use_float = use_float $
   , nonuniform = nonuniform 

; --- HANDLE EDGE CASE BETTER (AVOID WRAPS)
; --- AREA DECIMATION OF REGIONS DOESN'T HANDLE SHADOWING RIGHT NOW - FIX?

; --- ORDER OF GROWING (Z/XY)
; --- MAYBE ADD A VERBOSE FLAG AND GIVE SOME MORE OUTPUT?
; --- ADD SOME CAPABILITY TO CLIP HIGH NOISE REGIONS?

;+
; HISTORY:
;
;   20161212  LJ  add galaxy mask (inmask) as input
;   20170209  LJ  make min_pix/min_area to constrain lo_mask rather 
;                 than hi_mask
;   20170621  LJ  modify to calculate area of real 3D islands (i.e.  
;				  different regions in original 3D mask - "mask") rather 
;				  than the area of different regions of 2D mask - "mask_2d". 
;				  And also span the "hi_mask" to "lo_mask" first, then put on 
;                 the cariteria of "min_pix" and "min_area".
;   20170622  LJ  add a keyword "NONUNIFORM" to set the analysis based
;				  on "data/rms". If this keyword is not set, then the analysis is 
;				  default to based on real "data". However, whether the
;				  analysis is uniform or not depends on the input rms, so
;                 the rms cloud be uniform even when "NONUNIFORM" is set.
;   20171208  LJ  add output "outisland" - the assign cube for eash "island"
;   20180101  LJ  add mom0 and mom1 output "outmomzero" and "outmomone"
;
;-


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(infile) gt 0 then begin
     file_data = file_search(infile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Data not found.", /info
        return
     endif else begin
        cube = readfits(file_data, hdr)
     endelse
  endif

  if keyword_set(nonuniform) and (n_elements(rmsfile) eq 0) then begin
        message, "Rms file not found when NONUNIFORM is set", /info
        return
  endif


  if n_elements(rmsfile) gt 0 then begin
     file_rms = file_search(rmsfile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Noise not found.", /info
        return
     endif else begin
        rms = readfits(file_rms, rms_hdr)
     endelse
  endif

  if n_elements(indist) gt 0 then begin
     file_dist = file_search(indist, count=file_ct)
     if file_ct eq 0 then begin
        message, "Noise not found.", /info
        return
     endif else begin
        dist_cube = readfits(file_dist, dist_hdr)
     endelse
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET SOME DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(rms) eq 0 then $
     rms = mad(cube, /finite)

  if n_elements(hi_thresh) eq 0 then $
     hi_thresh = 5.0

  if n_elements(hi_nchan) eq 0 then $
     hi_nchan = 2

  if n_elements(lo_thresh) eq 0 then $
     lo_thresh = hi_thresh

  if n_elements(lo_nchan) eq 0 then $
     lo_nchan = hi_nchan

  if n_elements(clip_rms) eq 0 then $
     clip_rms = !values.f_nan

  sz = size(cube)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT THE DATA TO A SIGNIFICANCE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(verbose) then begin
     message, "Scaling cube by the noise.", /info
  endif

  rms_sz = size(rms)

; CLIP HIGH VALUES IN THE NOISE MAP
  if finite(clip_rms) then begin
     bad_rms = where(rms gt clip_rms, bad_ct)
     if bad_ct gt 0 then $
        rms[bad_rms] = !values.f_nan
  endif

; ... WE HAVE ONE NUMBER FOR THE RMS
  if keyword_set(nonuniform) then begin
	  if rms_sz[0] eq 0 then $
	     sig_cube = cube / rms

; ... WE HAVE AN RMS VECTOR (ASSUME Z)
	  if rms_sz[0] eq 1 then begin
	     sig_cube = cube
	     for i = 0, sz[3]-1 do $
	        sig_cube[*,*,i] = cube[*,*,i] / rms[i]
	  endif

; ... WE HAVE AN RMS MAP (ASSUME X-Y)
	  if rms_sz[0] eq 2 then begin
	     sig_cube = cube
	     for i = 0, sz[3]-1 do $
	        sig_cube[*,*,i] = cube[*,*,i] / rms
	  endif
	
; ... WE HAVE AN RMS CUBE
	  if rms_sz[0] eq 3 then begin
	     sig_cube = cube / rms
	  endif
  endif else begin
	  sig_cube = cube
  endelse

; ... IF DESIRED, INVERT THE CUBE
  if keyword_set(invert) then begin
     sig_cube *= -1.0
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE HIGH SIGNIFICANCE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(verbose) then begin
     message, "Making high threshold masks.", /info
  endif
  
; IDENTIFY ALL REGIONS WHERE nchan CHANNELS ARE ABOVE sig SIGMA
  conj = sig_cube gt hi_thresh
  for i = 1, hi_nchan-1 do $
     conj *= shift(sig_cube gt hi_thresh,0,0,i)

; SET ALL OF THE PIXELS IN THESE REGIONS TO 1 IN THE MASK
  for i = 1, hi_nchan-1 do $
     conj += shift(conj, 0,0,-1*i)
  hi_mask = conj ge 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE LOW SIGNIFICANCE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(verbose) then begin
     message, "Making low threshold masks.", /info
  endif
  
; IDENTIFY ALL REGIONS WHERE nchan CHANNELS ARE ABOVE sig SIGMA
  conj = sig_cube gt lo_thresh
  for i = 1, lo_nchan-1 do $
     conj *= shift(sig_cube gt lo_thresh,0,0,i)

; SET ALL OF THE PIXELS IN THESE REGIONS TO 1 IN THE MASK
  for i = 1, lo_nchan-1 do $
     conj += shift(conj, 0,0,-1*i)
  lo_mask = conj ge 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EXPAND THE HIGH SIGNIFICANCE CORE INTO THE LOW SIGNIFICANCE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(verbose) then begin
     message, "Expanding mask to lower threshold.", /info
  endif
  mask = grow_mask(hi_mask, constraint=lo_mask)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARE CONTIGUOUS REGIONS BY PIXELS OR AREA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(verbose) then begin
     message, "Paring low-volume regions from high threshold mask.", /info
  endif

; GET RID OF REGIONS SMALLER THAN A USER-SPECIFIED SIZE
  if n_elements(min_pix) gt 0 then begin
     reg = label_region(mask,/ulong)
     ind = where(reg ne 0, ct)
     if ct gt 0 then begin
        reg = reg[ind]
        max_reg = max(reg)
        for i = 1L, max_reg do begin
           if keyword_set(verbose) then $
              counter, i, max_reg, "Checking region "
           if (total(reg eq i) lt min_pix) then $
              mask[ind[where(reg eq i)]] = 0B     
        endfor
     endif
  endif

  ;;;;;;;; Modified by LJ ;;;;;;;;;;;;;;;
  if n_elements(min_edge_area) eq 0 then begin
    if n_elements(min_area) gt 0 then begin
       reg = label_region(mask, /ulong)
       ind = where(reg ne 0, ct)
       if ct gt 0 then begin
          max_reg = max(reg)
          sz = size(mask)
          for i = 1, max_reg do begin
             if keyword_set(verbose) then $
                counter, i, max_reg, "Checking region "
             index = where(reg eq i)
             xcor = index mod sz[1]
             ycor = (index mod (sz[1]*sz[2]))/sz[1]
             xyind = xcor + ycor * sz(1)
             this_area = n_elements(uniq(xyind, sort(xyind)))
             if (this_area lt min_area) then $
                mask[index] = 0B
          endfor
       endif
    endif  ; end min_area
  endif else begin
    if n_elements(min_area) gt 0 then begin
       reg = label_region(mask, /ulong)
       ind = where(reg ne 0, ct)
       if ct gt 0 then begin
          max_reg = max(reg)
          sz = size(mask)
          for i = 1, max_reg do begin
             if keyword_set(verbose) then $
                counter, i, max_reg, "Checking region "
             index = where(reg eq i)
             xcor = index mod sz[1]
             ycor = (index mod (sz[1]*sz[2]))/sz[1]
             xyind = xcor + ycor * sz(1)
             this_area = n_elements(uniq(xyind, sort(xyind)))
			 this_dist = min(dist_cube[index])

			 if this_dist le edge_pc then begin
               if (this_area lt min_area) then $
                 mask[index] = 0B
			 endif else begin
               if (this_area lt min_edge_area) then $
                 mask[index] = 0B
			 endelse

          endfor
       endif
    endif  ; end min_area
  endelse   ; end min_edge_area

  ;if n_elements(min_area) gt 0 then begin
  ;   mask_2d = total(lo_mask,3) gt 0
  ;   reg = label_region(mask_2d)
  ;   ind = where(reg ne 0, ct)
  ;   if ct gt 0 then begin
  ;      reg = reg[ind]
  ;      max_reg = max(reg)
  ;      for i = 1, max_reg do begin
  ;         if keyword_set(verbose) then $
  ;            counter, i, max_reg, "Checking region "        
  ;         if (total(reg eq i) lt min_area) then $
  ;            mask_2d[ind[where(reg eq i)]] = 0B     
  ;      endfor
  ;      for i = 0, (size(cube))[3]-1 do $
  ;         lo_mask[*,*,i] = lo_mask[*,*,i]*mask_2d
  ;   endif
  ;endif
  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;  GROW THE FINAL MASK IN THE XY OR Z DIRECTIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
  if n_elements(grow_xy) gt 0 then begin
     for i = 0, sz[3]-1 do $
        mask[*,*,i] =  grow_mask(mask[*,*,i], rad=grow_xy, /quiet)
  endif

  if n_elements(grow_z) gt 0 then begin
     for i = 0, grow_z do $
        mask = mask or shift(mask,0,0,i) or shift(mask,0,0,-i)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY A PRIOR, IF ONE IS SUPPLIED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(prior) gt 0 then begin
     mask *= prior
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(use_float) then $
     mask = float(mask)


  if n_elements(inmask) gt 0 then begin
     file_mask = file_search(inmask, count=file_ct)
     if file_ct eq 0 then begin
        message, "Galaxy Mask not found.", /info
        return
     endif else begin
		galaxy_mask=readfits(inmask, mask_hdr)
		inds=where(galaxy_mask gt 0, mask_ct)
		if mask_ct eq 0 then begin
			message, "Galaxy Mask is empty.",/info
		endif else begin	
			mask_copy=mask
			mask=mask*0.0
			mask[inds]=mask_copy[inds]
		endelse	
     endelse
  endif

  if n_elements(outfile) gt 0 then begin
	 hdr1 = hdr
     sxaddpar, hdr1, 'BUNIT', 'MASK'
     writefits, outfile, mask, hdr1
  endif

  if n_elements(outisland) gt 0 then begin
	 hdr2 = hdr
     sxaddpar, hdr2, 'BUNIT', 'MASK'
	 island_assign = label_region(mask,/ulong)
     writefits,outisland, island_assign, hdr2
  endif
  
; Can also output the mom0 fits
  cdelt3 = sxpar(hdr, 'CDELT3')
  mom0 = total(mask*cube, 3, /nan) * cdelt3
  new_bunit = sxpar(hdr, 'BUNIT') + ' km/s'
  print,'====================================================='
  print,'Total Flux in '+ new_bunit, total(mom0)
  print,'====================================================='
  if n_elements(outmomzero) gt 0 then begin
     cdelt3 = sxpar(hdr, 'CDELT3')
     mom0 = total(mask*cube, 3, /nan) * cdelt3
     new_bunit = sxpar(hdr, 'BUNIT') + ' km/s'
     hdr3 = hdr
     sxdelpar, hdr3, 'CTYPE3'
     sxdelpar, hdr3, 'CRVAL3'
     sxdelpar, hdr3, 'CRPIX3'
     sxdelpar, hdr3, 'CDELT3'
     sxdelpar, hdr3, 'CUNIT3'
     sxaddpar, hdr3, 'BUNIT', new_bunit
     sxaddpar, hdr3, 'DATAMAX', max(mom0,/nan)
     sxaddpar, hdr3, 'DATAMIN', min(mom0,/nan)
     writefits, outmomzero, mom0, hdr3
  endif

; Can also output the mom1 fits
  if n_elements(outmomone) gt 0 then begin
	 mom0 = total(mask*cube, 3, /nan)
	 sz = size(cube)
	 naxis3 = sxpar(hdr, 'NAXIS3')
	 crpix3 = sxpar(hdr, 'CRPIX3')
	 cdelt3 = sxpar(hdr, 'CDELT3')
	 crval3 = sxpar(hdr, 'CRVAL3')
	 cunit3 = sxpar(hdr, 'CUNIT3')
	 vaxis = (findgen(naxis3)+1 - crpix3)*cdelt3 + crval3
	 vcube = rebin(reform(vaxis,1,1,sz[3]), sz[1:3])
	 vcube = vcube * mask
	 mom1 = total(vcube*cube,3,/nan) / mom0
	 hdr4 = hdr
	 sxaddpar, hdr4, 'BUNIT', cunit3
     sxdelpar, hdr4, 'CTYPE3'
     sxdelpar, hdr4, 'CRVAL3'
     sxdelpar, hdr4, 'CRPIX3'
     sxdelpar, hdr4, 'CDELT3'
     sxdelpar, hdr4, 'CUNIT3'
     sxaddpar, hdr4, 'DATAMAX', max(mom1,/nan)
     sxaddpar, hdr4, 'DATAMIN', min(mom1,/nan)
     writefits, outmomone, mom1, hdr4
	  

  endif

end                             ; of make_cprops_mask
