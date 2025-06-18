Pro assign_check_convexity $
	, infile = infile $
	, inmask = inmask $
	, inassign = inassign $
	, rmsfile = rmsfile $
	, kernfile = kernfile $
	, minarea = minarea $
	, minpix = minpix $
	, minshape = minshape $
	, minconvexity = minconvexity $
	, mintexture = mintexture $
	, minUpar = minUpar $
	, delta = delta $
	, snr = delta_is_snr $
	, nonuniform = nonuniform $
	, outfile = outfile 


; NB - check whether every resolved cloud has convexity>=minconvexity and 
;      shape>=minshape. If not, then gradually move to smaller region (the
;      region constrained by upper contour_level), until the clouds has 
;      convexity>=minconvexity and shape>=minshape


; HISTORY
; Lijie Liu 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  if n_elements(infile) gt 0 then begin
     file_data = file_search(infile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Data not found.", /info
        return
     endif else begin
        data = readfits(file_data, hdr)
     endelse
  endif

  if n_elements(inmask) gt 0 then begin
     file_mask = file_search(inmask, count=file_ct)
     if file_ct eq 0 then begin
        message, "Mask not found.", /info
        return
     endif else begin
        mask = readfits(file_mask, hdr_m)
     endelse
  endif

  if n_elements(inassign) gt 0 then begin
     file_assign = file_search(inassign, count=file_ct)
     if file_ct eq 0 then begin
        message, "Inassign not found.", /info
        return
     endif else begin
        assign = readfits(file_assign, hdr_a)
     endelse
  endif


  if n_elements(delta_is_snr) eq 0 then delta_is_snr = 0

  if keyword_set(nonuniform) then begin
     if (n_elements(rmsfile) eq 0) then begin
        message, "Rms file not found when NONUNIFORM is set", /info
        return
     endif else begin
        file_rms = file_search(rmsfile, count=file_ct)
        if file_ct eq 0 then begin
            message, "Noise not found.", /info
            return
        endif else begin
           rms = readfits(file_rms, rms_hdr)
           data = data/rms
           sigma = 1.0

	       if delta_is_snr eq 0 then begin
		     delta = delta/rms
	       endif
        endelse
     endelse
  endif

  if not keyword_set(nonuniform) then begin
	  if delta_is_snr eq 1 then begin
         rms = readfits(file_rms, rms_hdr)
		 delta = delta * rms
	  endif
  endif

  ; CONVERT THE DELTA CUTOFF INTO A REAL INTENSITY CUTOFF (NECESSARY ATM)
  if delta_is_snr then $
     cutoff = sigma*delta $
  else $
     cutoff = delta

; PIXPERBEAM
  degperpix = sqrt(abs(sxpar(hdr, "cdelt1"))*abs(sxpar(hdr, "cdelt1")))
  bmaj_deg = sxpar(hdr, "BMAJ")
  bmin_deg = sxpar(hdr, "BMIN")
  beamfwhm_deg = sqrt(bmaj_deg*bmin_deg)
  beamfwhm_pix = beamfwhm_deg/degperpix
  pixperbeam = (beamfwhm_pix/2.0)^2*!pi/alog(2.)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE KERNELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; NOTE WHETHER WE GET THE MERGER MATRIX ALONG WITH THE KERNELS
  have_merger_matrix = 0B

  if n_elements(kernfile) gt 0 then begin
     readcol, kernfile, comment="#" $
              , format="L,L,L,F,F,F,F" $
              , kern_xpix, kern_ypix, kern_zpix $
              , kern_ra, kern_dec, kern_vel, kern_int
     xyv_to_ind, x=kern_xpix, y=kern_ypix, v=kern_zpix $
                 , sz=size(data), ind=kernel_ind
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFINITIONS AND DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; GET AN RMS ESTIMATE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if n_elements(sigma) eq 0 then begin
     sigma = mad(data,/finite)
     if keyword_set(nonuniform) then sigma = 1.0
  endif


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; GMC ASSIGNMENT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SET THE MINIMUM NUMBER OF LEVELS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if n_elements(nmin) eq 0 then $
     nmin = 100

; CALCULATE LEVELS TO WORK WITH
  vectorify, data $
             , mask = mask $
             , x = x, y = y, v = v, t = t $
             , ind = cubeindex

  szdata = size(data)
  cubify, x=x, y=y, v=v, t=t $
          , cube = cube $
          , pad = 3 $
          , dim_2d = (szdata[0] eq 2) $
          , indvec = cubeindex $
          , indcube = indcube

  spacing=0.05*sigma
  print, 'spacing in checking convexity', spacing ; added by Eric

  levels = $
     contour_values( $
     cube $
     , /linspace $
     , spacing=spacing $
     , nmin=nmin)

; GENERATE A CLOUD LIST IF NOT SUPPLIED
  ind = where(assign ne 0, ct)
  if ct eq 0 then begin
     message, "No valid assignments.", /info
     return
  endif

  if n_elements(cloudlist) eq 0 then begin
     cloudlist = assign[ind]
     cloudlist = cloudlist[sort(cloudlist)]
     cloudlist = cloudlist[uniq(cloudlist)]
  endif
  nclouds = n_elements(cloudlist)

; WORK OUT MERGER MATRIX
  if have_merger_matrix eq 0B then begin
     merger_matrix =  $
        mergefind_approx(data $
                         , kernel_ind $
                         , levels=levels $
                         , all_neighbors = all_neighbors $
                         , verbose = verbose)
  endif
  merger_matrix_nodiag = merger_matrix
  sz_merge = size(merger_matrix)
  ind = lindgen(sz_merge[1])
  merger_matrix_nodiag[ind,ind] = !values.f_nan

; UNIQUE LEVEL FOR EACH CLOUD
  unique_lev_array = make_array(nclouds) *!values.f_nan
  for i = 0, nclouds-1 do begin
	 if total(finite(merger_matrix_nodiag[i,*])) gt 0 then begin
        top_merger = max(merger_matrix_nodiag[i,*],/nan)
        unique_lev = min(levels[where(levels gt top_merger)],/nan)
     endif else begin
        unique_lev = min(levels,/nan)
     endelse

	 unique_lev_array[i] = unique_lev
  endfor

; VECTORIZE (SPEEDS UP SPARSE CASE)
  ind = where(assign ne 0, ct)
  assign_vec = assign[ind]
  t = data[ind]
  ind_to_xyv, ind, x=x, y=y, v=v, sz=size(data)


; LOOP OVER CLOUDS
  assign_new = assign * 0.
  cloudnr = 0
  for i = 0, nclouds-1 do begin
     if keyword_set(verbose) then begin
        counter, i+1, nclouds, "Moments for cloud "
     endif

	 this_ind = where(assign_vec eq cloudlist[i], ct)
     if ct eq 0 then continue

	 ; extract cloud pixel
     this_t = t[this_ind]
     this_x = x[this_ind]
     this_y = y[this_ind]
     this_v = v[this_ind]

	 ; construct minicube
	 pad = 2
     cubify $
     , x=this_x, y=this_y, v=this_v, t=this_t $
     , cube = minicube $
     , pad = pad $
     , location = location $
     , mask = cloudmask


     minisz = size(minicube)
	 xaxis = (indgen(minisz[1])+min(this_x)-pad)
     yaxis = (indgen(minisz[2])+min(this_y)-pad)
     xmap = (intarr(minisz[2])+1) ## xaxis
     ymap = yaxis ## (intarr(minisz[1])+1)


	 ; calculate texture
	 t0 = mean(this_t, /nan)
	 texture = sqrt(total((this_t - t0)^2.)/(n_elements(this_t)-1.))

	 ;    CALCULATE THE UPAR
	 mask1 = mask * 0.0
     mask1 = shift(cloudmask, 1,0,0)
     mask2 = shift(cloudmask, -1,0,0)
     mask3 = shift(cloudmask, 0,1,0)
     mask4 = shift(cloudmask, 0,-1,0)
     mask5 = shift(cloudmask, 0,0,1)
     mask6 = shift(cloudmask, 0,0,-1)
     mask_sum = mask1+ mask2 + mask3 + mask4 + $
       mask5 + mask6 + cloudmask


     inside = where(mask_sum eq 7, inside_ct)
     Vpar = n_elements(location)
     Qpar = Vpar - inside_ct ; number of the boundary pixels
     Upar = Qpar^3./(64.*!pi*Vpar^2.)


	 ; calculate convexity
     ; xmap and ymap
	 cloudmask2d = total(cloudmask,3,/nan)/(total(cloudmask, 3, /nan) > 1)
	 locat2d = where(cloudmask2d gt 0, area)   ; in pixel
	 vpix = n_elements(this_t)

	 x2d = xmap[locat2d]
	 y2d = ymap[locat2d]
	 convex_area = convexhull_area(x2d,y2d)
	 convexity2d = area/convex_area

	 mom0 = total(minicube,3,/nan)
	 I2d = mom0[locat2d]
     x2d = xmap[locat2d]
     y2d = ymap[locat2d]
     ;I2d = re_order(I2d)
     ;iVolume = total(I2d)  ; in [intensity_unit*pixel*pixel]
     ;convex_volume = convexhull_volume(x2d, y2d, I2d)
     ;convexity3d = iVolume/convex_volume
	 convexity3d = calc_convexity(x2d,y2d,I2d)

	 ;print,'cloudlist[i],shape,convexity3d,texture, Upar',cloudlist[i],convexity2d,convexity3d,texture, Upar

	 contrast = max(this_t) - min(this_t)

	 ; start loop
	 unique_lev = unique_lev_array[i]
	 bound_lev = min(this_t)
	

	 ; save original data cube
	 minicube0 = minicube
	 location0 = location
	 this_ind0 = this_ind
    ; print, i
	 while (area ge minarea) and (vpix ge minpix) and ((convexity2d lt minshape) or $
	 	 (convexity3d lt minconvexity) or (texture lt mintexture) or $
		 (Upar lt minUpar) or (contrast lt cutoff)) do begin

       bound_lev = bound_lev + spacing

	   ;this_ind = where((assign_vec eq cloudlist[i]) and (t ge bound_lev), ct)
	   reg = label_region(minicube0 ge bound_lev,/all_neighbors,/ULONG)
       max_ind = where(minicube0 eq max(minicube0,/nan))
       cube_ind = where(reg[location0] eq reg[max_ind[0]],this_ct)
       this_ind = this_ind0[cube_ind]
	   ;print,'max(reg),ct,this_ct',max(reg),ct,this_ct

	   ; extract cloud pixel
       this_t = t[this_ind]
       this_x = x[this_ind]
       this_y = y[this_ind]
       this_v = v[this_ind]

  	   ; construct minicube
  	   pad = 2
       cubify $
       , x=this_x, y=this_y, v=this_v, t=this_t $
       , cube = minicube $
       , pad = pad $
       , location = location $
       , mask = cloudmask
  
       ; xmap and ymap
       minisz = size(minicube)
  	   xaxis = (indgen(minisz[1])+min(this_x)-pad)
       yaxis = (indgen(minisz[2])+min(this_y)-pad)
       xmap = (intarr(minisz[2])+1) ## xaxis
       ymap = yaxis ## (intarr(minisz[1])+1)
  
       cloudmask2d = total(cloudmask,3,/nan)/(total(cloudmask, 3, /nan) > 1)
  	   locat2d = where(cloudmask2d gt 0, area)   ; in pixel

	   ; vpix
	   vpix = n_elements(this_t)
  

	   ; calculate texture
	   t0 = mean(this_t, /nan)
	   texture = sqrt(total((this_t - t0)^2.)/(n_elements(this_t)-1.))
  
       ;    CALCULATE THE UPAR
       mask1 = shift(cloudmask, 1,0,0)
       mask2 = shift(cloudmask, -1,0,0)
       mask3 = shift(cloudmask, 0,1,0)
       mask4 = shift(cloudmask, 0,-1,0)
       mask5 = shift(cloudmask, 0,0,1)
       mask6 = shift(cloudmask, 0,0,-1)
       mask_sum = mask1+ mask2 + mask3 + mask4 + $
           mask5 + mask6 + cloudmask

       inside = where(mask_sum eq 7, inside_ct)
       Vpar = n_elements(location)
       Qpar = Vpar - inside_ct ; number of the boundary pixels
       Upar = Qpar^3./(64.*!pi*Vpar^2.)

  	   ; calculate convexity
  	   x2d = xmap[locat2d]
  	   y2d = ymap[locat2d]
  	   convex_area = convexhull_area(x2d,y2d)
  	   convexity2d = area/convex_area
  
  	   mom0 = total(minicube,3,/nan)
	   I2d = mom0[locat2d]
       x2d = xmap[locat2d]
       y2d = ymap[locat2d]
       ;I2d = re_order(I2d)
       ;iVolume = total(I2d)  ; in [intensity_unit*pixel*pixel]
       ;convex_volume = convexhull_volume(x2d, y2d, I2d)
       ;convexity3d = iVolume/convex_volume
	   convexity3d = calc_convexity(x2d,y2d,I2d)

	   contrast = max(this_t) - min(this_t)

	 endwhile

    ; print, (convexity3d ge minconvexity)

	 ; check whether everything is okay
	 ;if area ge pixperbeam then $
	 ;  print,'id',cloudlist[i]
	 ;  print,'area,vpix,convexity2d,convexity3d,texture,Upar,contrast',$
	 ;     area,vpix,convexity2d,convexity3d,texture,Upar,contrast
	 ;  print, (area ge minarea) and (vpix ge minpix) and ((convexity2d lt minshape) or $
	 ;    (convexity3d lt minconvexity) or (texture lt mintexture) or $
	 ;	  (Upar lt minUpar)) and (contrast ge cutoff)


	 ; update assign_new
	 ;if area lt minarea then begin
	 ;  if (vpix ge minpix) and (contrast ge cutoff) then begin
	 ;    ;assign_new[ind[this_ind]] = cloudlist[i]
	 ;	 cloudnr = cloudnr+1
	 ;    assign_new[ind[this_ind]] = cloudnr
	 ;  endif
	 ;endif else begin
	 ;  if (vpix ge minpix) and (convexity2d ge minshape) $
	 ;	   and (convexity3d ge minconvexity) and $
	 ;	  (texture ge mintexture) and (Upar ge minUpar) and $
     ;    (contrast ge cutoff) then begin
	 ;      ;assign_new[ind[this_ind]] = cloudlist[i]
	 ;	   cloudnr = cloudnr+1
	 ;      assign_new[ind[this_ind]] = cloudnr
	 ;  endif
	 ;endelse
	 ; Keep All the CLouds
	 assign_new[ind[this_ind]] = cloudlist[i]

  endfor

; WRITE OUT OUTASSIGN
  if n_elements(hdr) gt 0 then begin
     out_hdr = hdr
     sxaddpar, out_hdr, "BUNIT", "ASSIGNMENT"
     writefits, outfile, assign_new, out_hdr
  endif else begin
     message, "Consider passing a header so that the assignment has astrometry.", /info
     writefits, outfile, assign_new
  endelse


	



End
