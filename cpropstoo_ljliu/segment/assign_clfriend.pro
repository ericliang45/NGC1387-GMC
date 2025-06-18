pro assign_clfriend $
   , kernels = kernel_ind $
   , kernfile = kernfile $
   , idlformat = idlformat $
   , data=data $
   , infile=infile $
   , mask=mask_in $
   , inmask=mask_file $
   , rmsfile=rmsfile $
   , spacing=spacing $
   , minlev=minlev $
   , nonuniform=nonuniform $
   , outfile=outfile $
   , diskmode=diskmode $
   , gal_par=gal_par

; NB - assumes equal weighting for velocity and spatial pixel; we may
;      want to add the ability to weight distance in spatial pixels
;      relative to distance in velocity pixels. In the case of a
;      weirdly sampled cube this is a big deal...

; Rosolowsky Clumpfind from CPROPSDIST Added by LIJIE LIU
;
; Added seeded/unseeded capability.
;  
; HISTORY
;
;

  compile_opt idl2

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

  if n_elements(mask_in) eq 0 then begin
     if n_elements(mask_file) gt 0 then begin
        full_file_mask = file_search(mask_file, count=file_ct)
        if file_ct eq 0 then begin
           message, "Mask file not found.", /info
           return
        endif else begin
           mask = readfits(full_file_mask, mask_hdr)
        endelse
     endif else begin
        message, "Defaulting to a mask of finite elements.", /info
        mask = finite(data)
     endelse
  endif else begin
     mask = mask_in
  endelse

  if keyword_set(diskmode) then begin
	  if n_elements(gal_par) ne 2 then begin
		  print,'ERROR: SET GAL_PAR FOR DISKMODE'
	  endif else begin
		  diskmode = 1
	      inc = gal_par[0]
	      pa = gal_par[1]
      endelse
  endif

; Return in the case of an empty mask
  if total(mask) eq 0 then begin
     message, "Empty mask. Returning.", /info
     return
  endif

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
        endelse
     endelse
  endif



; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE KERNELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(kernfile) gt 0 then begin
     if keyword_set(idlformat) then begin
        restore, kernfile
     endif else begin
        readcol, kernfile, comment="#" $
                 , format="L,L,L,F,F,F,F" $
                 , kern_xpix, kern_ypix, kern_zpix $
                 , kern_ra, kern_dec, kern_vel, kern_int
        xyv_to_ind, x=kern_xpix, y=kern_ypix, v=kern_zpix $
                    , sz=size(data), ind=kernel_ind
     endelse
  endif

  if n_elements(kernel_ind) ne 0 then begin
     seeded = 1B
  endif else begin
     seeded = 0B
  endelse

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFINITIONS AND DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; IF NO SPACING IS DEFINED DEFAULT TO 2 TIMES THE RMS NOISE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if n_elements(sigma) eq 0 then begin
     sigma = mad(data,/finite)
     if keyword_set(nonuniform) then sigma = 1.0
  endif


  if n_elements(spacing) eq 0 then begin
     spacing = 0.2*sigma  
  endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SET THE MINIMUM LEVEL TO CONSIDER
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if n_elements(minlev) eq 0 then $
     minlev = min(data[where(mask)], /nan)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARE CUBE TO MINIMUM SIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Convert the data inside the mask to a vector

  vectorify, data $
             , mask = mask $
             , x = x, y = y, v = v, t = t $
             , ind = cubeindex

; Rebuild a minimum-sized cube from the vectorized data.

  szdata = size(data)
  cubify, x=x, y=y, v=v, t=t $
          , cube = minicube $
          , pad = 3 $
          , dim_2d = (szdata[0] eq 2) $
          , indvec = cubeindex $
          , indcube = indcube
  
  if seeded then begin
     minikern = kernel_ind
     for i = 0, n_elements(minikern)-1 do $
        minikern[i] = where(indcube eq kernel_ind[i])
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ASSIGNMENT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; Rosolowsky Clumpfind
; Calculate ppbeam
  degperpix = sqrt(abs(sxpar(hdr, "cdelt1"))*abs(sxpar(hdr, "cdelt1")))
  bmaj_deg = sxpar(hdr, "BMAJ")
  bmin_deg = sxpar(hdr, "BMIN")
  beamfwhm_deg = sqrt(bmaj_deg*bmin_deg)
  beamfwhm_pix = beamfwhm_deg/degperpix
  pixperbeam = (beamfwhm_pix/2.0)^2*!pi/alog(2)


; Assign clump/GMC using Rosolowsky Clumpfind way
; Note the output assignment (astr.assignment) count from 0 for the first
; kernel
  cube_partition, minicube $
	  , level = spacing $
	  , terminate = 0.0, astr = astr $
	  , kernels = minikern, ppbeam = pixperbeam $
	  , diskmode = diskmode, gal_par = gal_par

  sz = size(minicube)
  miniassign = lonarr(sz[1], sz[2], sz[3])
  asgn_ind = astr.index
  miniassign[asgn_ind] = astr.assignment+1
  

; INITIALIZE ASSIGNMENT CUBE
  sz = size(data)
  assign = lonarr(sz[1], sz[2], sz[3])
  assign_ind = where(miniassign gt 0)
  assign[indcube[assign_ind]] = miniassign[assign_ind]


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; OUTPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if n_elements(hdr) gt 0 then begin
     out_hdr = hdr
     sxaddpar, out_hdr, "BUNIT", "ASSIGNMENT"
     writefits, outfile, assign, out_hdr
  endif else begin
     message, "Consider passing a header so that the assignment has astrometry.", /info
     writefits, outfile, assign
  endelse

end


end

