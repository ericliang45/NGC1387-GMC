function alllocmax $
   , cubein $
   , indcube = indcube $
   , friends = friends $
   , specfriends = specfriends $
   , search_kern = search_kern $
   , exact = exact $
   , patchmode = patchmode $
   , patch_friends = patch_friends $
   , patch_specfriends = patch_specfriends $
   , verbose = verbose

;+
;
; NAME:
;
;   ALLLOCMAX()
;
; PURPOSE:
;   To establish a candidate set of local maxima within a data cube
;   by searching over a 2F+1 by 2F+1 by 2S+1 box of pixels for a point
;   that's larger than all others in the box.  F and S default to 1.
;
; CALLING SEQUENCE:
;   local_maxima = ALLLOCMAX(cube, [friends = friends, specfriends =
;   specfriends, indcube = indcube])
;
; INPUTS:
;   CUBE -- data cube to find local maxima in.
;             
; KEYWORD PARAMETERS:
;   INDCUBE -- (optional) if supplied, then the indices returned are
;              drawn from this cube instead of indexing the supplied
;              cube itself.
;   FRIENDS -- Sets the search region in the position dimensions to
;              be 2*FRIENDS+1 in size.  
;   SPECFRIENDS -- Sets the search region in the velocity dimensions
;              to be 2*SPECFRIENDS+1 in size.  
;
; OUTPUTS:
;   LOCAL_MAXIMA -- indices in CUBE (or in INDCUBE if set) that are
;                   local maxima.
;
; MODIFICATION HISTORY:
;
;       Documented -- Fri Sep 2 15:48:17 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;       Fixed bug concering NaNs in 'quick' dilate kernel search
;                     Oct 5 2014, Andreas Schruba <schruba@mpe.mpg.de>
; 
;  20180129     Lijie Liu    Add patch mode
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS AND DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(cubein)

  if sz[0] eq 2 then begin
     specfriends = 0
     sz[3] = 1
  endif

; INITIALIZE THE DEFAULT BOX SIZE  
  if n_elements(search_kern) eq 0 then begin

     if (n_elements(friends) eq 0) then $
        friends = 1

     if n_elements(specfriends) eq 0 then $
        specfriends = 1 

;    INITIALIZE A TRIVIAL SEARCH KERNEL
     search_kern = bytarr(friends*2+1, friends*2+1, specfriends*2+1)+1B     

  endif else begin

;    PROCESS A USER-SUPPLIED SEARCH KERNEL
     sz_search = size(search_kern)

     if (sz_search[1] mod 2 ne 1) or $
        (sz_search[2] mod 2 ne 1) then begin
        message, "Search kernel must have odd dimensions.", /info
        return, !values.f_nan
     endif
     friends = ((sz_search[1] > sz_search[2]) - 1) / 2
     
     if sz_search[0] eq 3 then begin
        if (sz_search[3] mod 2 ne 1) then begin
           message, "Search kernel must have odd dimensions.", /info
           return, !values.f_nan
        endif 
        spec_friends = (sz_search[3] - 1) / 2
     endif

  endelse

; INITIALIZE THE PATCH KEYWORDS 
  if keyword_set(patchmode) then begin
	  if n_elements(patch_friends) eq 0 then $
		  patch_friends = friends

	  if n_elements(patch_specfriends) eq 0 then $
		  patch_specfriends = specfriends
  endif

; DEFAULT QUANTIZATION
  nquant = 100000L

;
; A local maximum is defined to be a point greater than all the other
; points around it within an N x N x M box, where N and M are defined
; by the friends and specfriends keywords. We find these points in one
; of two ways. The quick way is to leverage IDL's native DILATE
; routine. The downside of this is that this requires quantizing the
; cube (although we do this at a very coarse level). The exact
; alternative uses SHIFT to roll the cube along each axis. This can be
; much slower for very large cubes.
;

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FAST: QUANTIZE THE CUBE AND USE IDL'S IMAGE DILATION ROUTINE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Identify a set of maxima by quantizing the cube, applying image
; dilation, and then comparing the dilated cube to the original.

  if keyword_set(exact) eq 0 then begin

;    RECAST THE CUBE AS UNSIGNED LONGS     
     maxval = max(cubein,/nan, min=minval)
     cube = ulong((cubein - minval)/(maxval - minval)*nquant)
     
;    MAY NOT BE NECESSARY... BUT AVOID NANS
     bad_ind = where(finite(cubein) eq 0, bad_ct)
     if bad_ct gt 0 then $
        cube[bad_ind] = 0

;    APPLY DILATION OPERATOR
     dilated_cube = $
        dilate(cube, search_kern, /gray, /preserve)
     
     lmaxcube = cube eq dilated_cube and dilated_cube ne 0

     print, 'NUMBER OF LMAX *before* PATCHMODE:', total(lmaxcube)

	 if keyword_set(patchmode) then begin
		 cubein_float = cubein
		 bad_ind = where(finite(cubein) eq 0, ct)
		 cubein_float[bad_ind] = 0.0
		 patch_cubein = cubein_float * 0.0
		if patch_specfriends gt 0 then begin
    	    for k = -patch_specfriends, patch_specfriends do begin
    	      for j = -patch_friends, patch_friends do begin
                for i = -patch_friends, patch_friends do begin
                  patch_cubein = patch_cubein + shift(cubein_float, i, j, k)
                endfor
              endfor
            endfor
		endif else begin
    	    for j = -patch_friends, patch_friends do begin
              for i = -patch_friends, patch_friends do begin
                patch_cubein = patch_cubein + shift(cubein_float, i, j, 0)
              endfor
            endfor
	    endelse

         maxval = max(patch_cubein, /nan, min = minval)
         patch_cube = ulong((patch_cubein - minval)/(maxval - minval)*nquant)
         
         bad_ind = where(finite(patch_cubein) eq 0, bad_ct)
         if bad_ct gt 0 then $
                patch_cube[bad_ind] = 0

         patch_search_kern = bytarr(3*(2*patch_friends+1), 3*(2*patch_friends+1), $
                    (2*patch_specfriends+1))+0B
         x_patch = [patch_friends, 3*patch_friends+1, 5*patch_friends+2]
         y_patch = [patch_friends, 3*patch_friends+1, 5*patch_friends+2]
         for i = 0, 2 do begin
           for j = 0, 2 do begin
             patch_search_kern[x_patch[i], y_patch[j], patch_specfriends] = 1B
           endfor
         endfor

		 ; We consider a (2*patch_friends+1) x (2*patch_friends+1) patch centered 
		 ; at the each pixel and compared this patch with all its eight neighboring 
		 ; (2*patch_friends+1) x (2*patch_friends+1) patches to determine 
		 ; whether it is a local maximum. The sums of the pixle values in each 
		 ; patch are compared.
	     dilated_patch_cube = $
    	    dilate(patch_cube, patch_search_kern, /gray, /preserve)
     
	     patch_lmaxcube = patch_cube eq dilated_patch_cube $
				and dilated_patch_cube ne 0

         bad_ind = where(finite(cubein) eq 0, bad_ct)
         if bad_ct gt 0 then $
                patch_lmaxcube[bad_ind] = 0

         lmaxcube = lmaxcube * patch_lmaxcube
         print, 'NUMBER OF LMAX *after* PATCHMODE:', total(lmaxcube)

	 endif  ; end keyword patchmode

  endif  ; end keyword exact

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EXACT: IDENTIFY LOCAL MAXIMA BY ROLLING THE CUBE ALONG EACH AXIS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Identify a set of local maxima by rolling ("shift"ing) the cube
; along each axis comparing these to the original cube.

  if keyword_set(exact) then begin

;    COPY CUBE
     cube = cubein
     
;    INITIALIZE EVERYTHING TO BE A LOCAL MAX, WE WILL REJECT THEM OVER TIME
     lmaxcube = bytarr(sz[1], sz[2], sz[3]) + 1B

;    SET NON-FINITE VALUES TO NEGATIVE INFINITY
     badind = where(finite(cubein) eq 0, badct)
     if (badct gt 0) then begin
        cube[badind] = -!values.f_infinity
        lmaxcube[badind] = 0B
     endif

;    INITIALIZE COUNTER FOR OUTPUT
     total_rolls = (2*specfriends+1)*(2*friends+1)^2
     roll_count = 0

;    ROLL THE CUBE
     if specfriends gt 0 then begin
        for k = -specfriends, specfriends do begin
           for j = -friends, friends do begin
              for i = -friends, friends do begin

;                IF REQUESTED, GIVE THE USER A PROGRESS METER
                 roll_count += 1
                 if keyword_set(verbose) then begin
                    counter, roll_count, total_rolls, "Cube shift "
                 endif

;                DON'T COMPARE TO SELF
                 if i eq 0 and j eq 0 and k eq 0 then $
                    continue

;                COMPARE ONLY IN AREA DEFINED BY SEARCH KERNEL
                 if search_kern[i+friends,j+friends,k+specfriends] eq 0B then $
                    continue
                 
;                CHECK IF ORIGINAL CUBE EXCEEDS THE SHIFTED CUBE
                 ;lmaxcube *= (cube gt shift(cube, i, j, k))
                 lmaxcube *= (cube ge shift(cube, i, j, k))   ; LJ 

              endfor
           endfor
        endfor
     endif else begin
        for j = -friends, friends do begin
           for i = -friends, friends do begin

;             IF REQUESTED, GIVE THE USER A PROGRESS METER
              roll_count += 1
              if keyword_set(verbose) then begin
                 counter, roll_count, total_rolls, "Cube shift "
              endif

;             DON'T COMPARE TO SELF
              if i eq 0 and j eq 0 then $
                 continue

;             COMPARE ONLY IN AREA DEFINED BY SEARCH KERNEL
              if search_kern[i+friends,j+friends] eq 0B then $
                 continue
              
;             CHECK IF ORIGINAL CUBE EXCEEDS THE SHIFTED CUBE
              lmaxcube *= (cube gt shift(cube, i, j))

           endfor
        endfor
     endelse
     print, 'NUMBER OF LMAX *before* PATCHMODE:', total(lmaxcube)

;    ROLL THE PATCH CUBE
	 if keyword_set(patchmode) then begin
		cubein_float = cubein
		bad_ind = where(finite(cubein) eq 0, ct)
		cubein_float[bad_ind] = 0.0
		patch_cubein = cubein_float * 0.0
		if patch_specfriends gt 0 then begin
    	    for k = -patch_specfriends, patch_specfriends do begin
    	      for j = -patch_friends, patch_friends do begin
                for i = -patch_friends, patch_friends do begin
                  patch_cubein = patch_cubein + shift(cubein_float, i, j, k)
                endfor
              endfor
            endfor
		endif else begin
    	    for j = -patch_friends, patch_friends do begin
              for i = -patch_friends, patch_friends do begin
                patch_cubein = patch_cubein + shift(cubein_float, i, j)
              endfor
            endfor
	    endelse

		patch_cube = patch_cubein
     
		patch_lmaxcube = bytarr(sz[1], sz[2], sz[3]) + 1B

		badind = where(finite(patch_cubein) eq 0, badct)
		if (badct gt 0) then begin
	        patch_cube[badind] = -!values.f_infinity
			patch_lmaxcube[badind] = 0B
	    endif

		; We consider a (2*patch_friends+1) x (2*patch_friends+1) patch centered at the
		; each pixel and compared this patch with all its eight (2*patch_friends+1)
		; x (2*patch_friends+1) neighboring patches to determine whether it is a
		; local maximum. The sums of the pixle values in each patch are
		; compared.
        for j = -1, 1 do begin
           for i = -1, 1 do begin

;             DON'T COMPARE TO SELF
              if i eq 0 and j eq 0 then $
                 continue

              
;             CHECK IF ORIGINAL CUBE EXCEEDS THE SHIFTED CUBE
              patch_lmaxcube *= (patch_cube ge shift(patch_cube, $
				  i*(2*patch_friends+1), j*(2*patch_friends+1), 0))

           endfor
        endfor

        bad_ind = where(finite(cubein) eq 0, bad_ct)
        if bad_ct gt 0 then $
               patch_lmaxcube[bad_ind] = 0

        lmaxcube = lmaxcube * patch_lmaxcube
        print, 'NUMBER OF LMAX *after* PATCHMODE:', total(lmaxcube)

	  endif  ; end keyword patchmode


  endif   ; end keyword exact
  
; EXTRACT INDICES OF LOCAL MAXIMA
  lmaxind = where(lmaxcube eq 1B, num)
  if (num eq 0) then begin
     message, 'No true local max found, defaulting to high point in data.', /con
     dummy = max(lmaxcube, lmaxind, /nan)
  endif

; IF THE INDEX CUBE IS SUPPLIED AND THERE ARE LOCAL MAXIMA THEN
; SUBSTITUTE THE INDICES FROM THE CUBE FOR THE ACTUAL INDICES
  if ((n_elements(indcube) gt 0)) then begin
     lmaxind = indcube[lmaxind]
  endif 

  return, lmaxind
end                             ; OF ALLLOCMAX


