pro find_local_max_loop $
   , data=data $
   , infile=infile $
   , mask=mask_in $
   , inmask=mask_file $
   , trueval=trueval $
   , hdr=hdr $
   , friends = friends $
   , specfriends = specfriends $
   , minpix = minpix $
   , minvchan = minvchan $
   , minarea = minarea $
   , minshape = minshape $
   , minconvexity = minconvexity $
   , mintexture = mintexture $
   , patchmode = patchmode $
   , fscale = fscale $
   , sigdiscont = sigdiscont $
   , nredun = nredun $
   , all_neighbors = all_neighbors $
   , kernels = sub_kernels $ ; modified by Eric ; kernels = kernels
   , delta = delta $
   , snr = delta_is_snr $
   , rmsfile = rmsfile $
   , inrms = rms $
   , minval = minval $
   , nodecimate = nodecimate $
   , no_deriv_decimate = no_deriv_decimate $
   , justdecimate = justdecimate $
   , idl_out = idl_out $
   , text_out = text_out $
   , verbose = verbose $
   , bclip = bclip $
   , nonuniform = nonuniform


;+
;
; NAME:
;
;   FIND_LOCAL_MAX_LOOP
;
; PURPOSE:
; Select kernels for big clouds to small clouds. Normally, the values of 
; minarea/minpix are set in a dereasing order. So it first selects all the 
; clouds have biggest values of minarea/minpix which satisfy all the selecting 
; requirements (minarea/minpix/minvchan/delta) and meanwhile do not have 
; substrucutres (minshape/minconvexity).  Once the big clouds are selected, 
; then only keep one original local maxima and decimate all the other local 
; maxima inside the cloud. The kernels outside of the big clouds are kept 
; together with selcted kernels for big clouds, and continue to be decimated 
; in next run.
;
; CALLING SEQUENCE:
;    
;
; INPUTS:
;   INFILE -- Path to a .fits cube.
;   DATA -- (optional) CO Cube.  
;   INMASK -- Path to .fits mask (must be same size as data)    
;   MASK -- (optional) CO byte mask.
;   TRUEVAL -- (optional) truth value of mask. Defaults to >= 1.
;   HDR -- (optional) .fits Header (required if no filepath is specified). 
;   text_out -- Text file to save kernel locations in ascii format.
;   idl_out -- IDL file (.sav or .idl) to save kernel
;              locations in IDL format (variable name = kernel_ind).
;             
; KEYWORD PARAMETERS:
;   FRIENDS -- (optional) Pixels to search over in the x-y plane. Total
;              search box length is 2*Friends+1. Default is Friends=3
;   SPECFRIENDS -- (optional) Pixels to search over in the v plane.
;              Total search box length is 2*Specfriends+1. Default is
;              Specfriends=1
;
; OUTPUTS: 
;   KERNELS -- Array of local maxima.

; MODIFICATION HISTORY:
;      Originally written by Adam Leroy and Erik Rosolowsky.
;
;      Some documentation -- Mon Nov 25, 2013  Stephen Pardy 
;                     <spardy@astro.wisc.edu>
; 
;-



;+
;
; TBD:
;
; - proper island handling in the decimation
; - default mask creation call? debateable
; - noise assumed homogeneous if rmsfile unassigned
;
;   20170214    LJ  add derivative decimation based on cprops_dist
;   20170622    LJ  add a keyword "NONUNIFORM" to set the analysis based
;                   on "data/rms". If this keyword is not set, then the analysis is 
;                   default to based on real "data". However, whether the
;                   analysis is uniform or not depends on the input rms, so
;                   the rms cloud be uniform even when "NONUNIFORM" is set.
;   20170830    LJ  add a keyword "bclip" to rescale the data to reduce the
;					effect of substructure, more details refer to Roslowsky
;					& Leroy (2006)
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

; If requested, use only the mask where it equals a certain true
; value. Useful for analyzing only part of an assignment cube, for
; example.
  if n_elements(trueval) ne 0 then begin
     mask = mask eq trueval
  endif

; Return in the case of an empty mask
  if total(mask) eq 0 then begin
     message, "Empty mask. Returning.", /info
     return
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS AND DEFINTIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  szdata = size(data)

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ERROR TRAP
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if n_elements(kernels) eq 0 and $
     keyword_set(just_decimate) then begin
     message, "Cannot decimate without a kernel list.", /info
  endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SIGNAL-TO-NOISE AND RMS UNITS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; DEFINE UNITS OF DELTA (S/N OR REAL)
  if n_elements(delta_is_snr) eq 0 then $
     delta_is_snr = 0B


  if keyword_set(nonuniform) and (n_elements(rmsfile) eq 0) then begin
        message, "Rms file not found when NONUNIFORM is set", /info
        return
  endif

  if keyword_set(nonuniform) and (delta_is_snr eq 0B) then begin
        message, "Keyword SNR need to be set when NONUNIFORM is set", /info
        message, "to make sure delta is given in unit of sigma !", /info
        return
  endif


  if n_elements(rmsfile) gt 0 then begin
     file_rms = file_search(rmsfile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Noise file not found.", /info
        return
     endif else begin
        rms = readfits(file_rms, rms_hdr)
      ; In this case RMS has to be a CUBE
        sz_rms  = size(rms)
        sz_data = size(data)
        if sz_rms[0] eq 2 and sz_data[0] eq 3 then $
           rms = rebin(reform(rms,sz_rms[1],sz_rms[2],1),sz_data[1:3])
     endelse
  endif else begin
     if n_elements(rms) eq 0 then begin        
        rms = mad(data, /finite)
     endif
  endelse

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SEARCH CRITERIA ...
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; SPATIAL SEARCH AREA IN PIXELS
  if n_elements(friends) eq 0 then $
     friends = 3

; THE SEARCH AREA ALONG THE THIRD DIMENSION IN PIXELS
  if n_elements(specfriends) eq 0 then $
     specfriends = 1
  if szdata[0] eq 2 then $
     specfriends = 0

; DEFINE CONNECTEDNESS
  if not(keyword_set(all_neighbors)) then $
     all_neighbors = 0b

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; REJECTION CRITERIA ...
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; CONTRAST WITH MERGE LEVEL
  if n_elements(delta) eq 0 then $
     delta = 3

; MINIMUM VOLUME
  if n_elements(minval) eq 0 then $
     minval = 0

; MINIMUM VOLUME
  if n_elements(minpix) eq 0 then $
     minpix = [4]

; MINIMUM AREA
  if n_elements(minarea) eq 0 then $
     minarea = [1]

  if n_elements(minarea) ne n_elements(minpix) then begin
	  print,'ERROR: THE SIZE OF MINAREA AND MINPIX ARE NOT EQUAL!'
	  return
  endif


; MINIMUM SPECTRAL EXTENT
  if n_elements(minvchan) eq 0 then $
     minvchan = 1

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DEFAULT OUTPUT NAMES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if n_elements(text_out) eq 0 then $
     text_out = "lmax.txt"

;  if n_elements(idl_out) eq 0 then $
;     idl_out = "lmax.idl"      ;do not generate lmax.idl by LJ

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARE CUBE TO MINIMUM SIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Convert the data inside the mask to a vector. Do the same for the
; noise if we have a cube.

; Set bclip to rescale the data to reduce the effect of substructure.
; More details can be found in Rosolowsky & Leroy (2006).
; However, using typical values Tclip = 2.5 K has little influence 
; on extragalactic data where beam deconvolution typically averages out 
; the presence of bright substructure within the clouds.
  if n_elements(bclip) gt 0 then begin
	  clip_index = where(data gt bclip, clip_count)
	  if clip_count gt 0 then begin
		  orig_data_clip = data[clip_index]
	      data[clip_index] = $
		  (bclip*(1+atan(data[clip_index]/bclip-1))) 
	  endif
  endif else clip_count = 0

 
  if keyword_set(nonuniform) then begin
	  data = data/rms
	  sigma = 1
  endif
      
 

  vectorify, data $
             , mask = mask $
             , x = x, y = y, v = v, t = t $
             , ind = cubeindex

  if n_elements(rms) eq n_elements(data) then begin
     vectorify, rms $
                , mask = mask $
                , x = x, y = y, v = v, t = e_t $
                , ind = cubeindex     
  endif else begin
     if n_elements(rms) ne 1 then begin
        message, "Size of rms is unexpected.", /info
        return
     endif
  endelse

; Rebuild a minimum-sized cube from the vectorized data.

  cubify, x=x, y=y, v=v, t=t $
          , cube = minicube $
;         , pad = (friends > specfriends) $
          , pad = 3 $
          , dim_2d = (szdata[0] eq 2) $
          , indvec = cubeindex $
          , indcube = indcube


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY KERNELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; This step runs unless the "justdecimate" flag is set, in which case
; it is assumed that kernels have been supplied by the user.

  if keyword_set(justdecimate) eq 0 then begin
     sub_kernels = alllocmax(minicube $
                             , friends = friends $
                             , specfriends = specfriends $
							 , patchmode = patchmode $
                             , verbose = verbose)
  endif

; This step runs unless the "nodecimate" flag is set. The subroutine
; rejects all but significant kernels that are defined by a set of
; user-tunable quantities: delta, minval, minpix, minarea, minvchan.

  if keyword_set(nodecimate) eq 0 then begin
     if not keyword_set(nonuniform) then begin
     	sigma = mad(data,/finite) ; estimate here, since RMS of masked cube may be different
	 endif else begin
		 sigma = 1.
	 endelse
;    print,"Calling decimate_kernels with sigma value: ", sigma
;    for first few minarea/minpix only get rid of kernels inside
;    each unique_level of valid kernels
	 j = 0
     for i = 0, n_elements(minarea)-2 do begin
       sub_kernels = $
          decimate_kernels_inside(sub_kernels $
                           , minicube $
                           , all_neighbors = all_neighbors $
                           , delta = delta $
                           , snr = delta_is_snr $
                           , sigma = sigma $
                           , minval = minval $
                           , minpix = minpix[i] $
                           , minarea = minarea[i] $
                           , minvchan = minvchan $
						   , minshape = minshape $
						   , minconvexity = minconvexity $
						   , mintexture = mintexture $
  						   , fscale = fscale $
  					       , sigdiscont = sigdiscont $
  						   , nredun = nredun $
  						   , no_deriv_decimate = no_deriv_decimate $
                           , verbose = verbose $
                           , valid_merger = merger_matrix)
		j = i+1
		message, 'KERNELS LEFT FOR THE RUN '+strtrim(i+1)+': '+strtrim(n_elements(sub_kernels)), /con
	 endfor				   
;    For last minarrea/minpix get rid of all the kernels not valid
;    Not to consider minshape/minconvexity for last run, because the
;    shapes of each kernel are not independent, they are associated with
;    each other. After assignment, we are then in the right place to 
;    remove the GMC with weired shape/convexity.
     sub_kernels = $
        decimate_kernels(sub_kernels $
                         , minicube $
                         , all_neighbors = all_neighbors $
                         , delta = delta $
                         , snr = delta_is_snr $
                         , sigma = sigma $
                         , minval = minval $
                         , minpix = minpix[j] $
                         , minarea = minarea[j] $
                         , minvchan = minvchan $
						 ;, minshape = minshape $
						 ;, minconvexity = minconvexity $
  	  				     , fscale = fscale $
  	  			         , sigdiscont = sigdiscont $
  	  				     , nredun = nredun $
  	  				     , no_deriv_decimate = no_deriv_decimate $
                         , verbose = verbose $
                         , valid_merger = merger_matrix)
	 message, 'KERNELS LEFT FOR THE RUN '+strtrim(i+1)+': '+strtrim(n_elements(sub_kernels)), /con

  endif

  kernels = indcube[sub_kernels]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  write_kernels $
     , kernels $
     , sz = size(data) $
     , cube = data $
     , hdr = hdr $
     , text_file = text_out $
     , idl_file = idl_out $
     , merger = merger_matrix
 
 
  return

end
