pro cube_to_moments $
   , data = data $
   , infile = infile $
   , assign = assign $
   , inassign = inassign $
   , rms = rms $
   , rmsfile = rmsfile $
   , hdr = hdr $
   , outfile = outfile $
   , verbose = verbose $
   , cloudlist = cloudlist $
   , bootstrap = bootstrap


;   20170213  LJ  add bootstrap uncertainties based on cprops_dist

  compile_opt idl2
   
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFINE THE SET OF MODULES THAT WE WILL WORK WITH
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; Each modules adds fields and calculations to the moment structure.

  modules = [ "classic" $
            , "gausscorr" $
            , "area" ]

  if n_elements(extra_modules) gt 0 then begin
     modules = [modules, extra_modules]
  endif
  n_mod = n_elements(modules)


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(infile) gt 0 then begin
     file_data = file_search(infile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Data not found.", /info
        return
     endif else begin
        data = readfits(file_data, hdr, /silent)
     endelse
  endif

  ; extract beam
  extast, hdr, astr
  xy2ad, [0,1], [0,0], astr, ra, dec
  degperpix = sphdist(ra[0], dec[0], ra[1], dec[1], /deg)
  bmaj_deg = sxpar(hdr, "BMAJ") 
  bmin_deg = sxpar(hdr, "BMIN") 
  bmaj_pix = bmaj_deg / degperpix
  bmin_pix = bmin_deg / degperpix
  bpa_deg = sxpar(hdr, "BPA")
  beamfwhm_pix = sqrt(bmaj_pix * bmin_pix)
  ppbeam = !pi * (beamfwhm_pix/2.)^2./alog(2.)
  indfac = sqrt(ppbeam)
;   print, indfac
; end

  if n_elements(assign) eq 0 then begin
     file_assign = file_search(inassign, count=file_ct)
     if file_ct eq 0 then begin
        message, "Assignment cube not found.", /info
        return
     endif else begin
        assign = readfits(file_assign, assign_hdr, /silent)
     endelse
  endif
  
  if n_elements(rms) eq 0 then begin
     file_ct  = 0
     if n_elements(rmsfile) gt 0 then $
       file_rms = file_search(rmsfile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Noise cube not found. Not critical.", /info
     endif else begin
        rms = readfits(file_rms, rms_hdr, /silent)
     endelse
  endif

  if n_elements(bootstrap) eq 0 then begin
	  bootstrap = 0
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALL THE MEASUREMENT CODE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; EXTRACT PIXELS WITH ASSIGNMENTS
  ind = where(assign ne 0, ct)
  if ct eq 0 then begin
     message, "No valid assignments.", /info
     return
  endif  

; GENERATE A CLOUD LIST IF NOT SUPPLIED
  if n_elements(cloudlist) eq 0 then begin
     cloudlist = assign[ind]
     cloudlist = cloudlist[sort(cloudlist)]
     cloudlist = cloudlist[uniq(cloudlist)]
  endif

; VECTORIZE (SPEEDS UP SPARSE CASE)
  assign_vec = assign[ind]
  t = data[ind]
  ind_to_xyv, ind, x=x, y=y, v=v, sz=size(data)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALL THE MEASUREMENT CODE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; LOOP OVER CLOUDS
  nclouds = n_elements(cloudlist)

  bootmomra=[]
  for i = 0, nclouds-1 do begin
       
     if keyword_set(verbose) then begin
        counter, i+1, nclouds, "Moments for cloud "
     endif

     ind = where(assign_vec eq cloudlist[i], ct)
     if ct eq 0 then continue    

     this_t = t[ind]
     this_x = x[ind]
     this_y = y[ind]
     this_v = v[ind]

     
     this_mom = $
        measure_moments( x=this_x, y=this_y, v=this_v, t=this_t $
                       , /extrap, extarg=0, hdr = hdr $
                       , empty_props = empty_props, rms = rms )

	 ; Calculate Uncertainties
	 if (bootstrap gt 0) then begin   
	   npts=n_elements(ind)

       mom2maj_extrap_deconv = findgen(bootstrap)
       mom2min_extrap_deconv = findgen(bootstrap)
       momposang_deconv = findgen(bootstrap)
       mom2v_extrap_deconv = findgen(bootstrap)
       mom2v_gauss_extrap_deconv = findgen(bootstrap)
       mom0_extrap = findgen(bootstrap)
	   for j = 0, bootstrap-1 do begin
		 ;GENERATE A NEW SET OF DATA POINTS FROM BOOTSTRP 
		 ;RESAMPLING FROM CLOUD DATA
		 bootind=fix(randomu(seed, npts)*npts)
		 bootx=this_x[bootind]
		 booty=this_y[bootind]
		 bootv=this_v[bootind]
		 boott=this_t[bootind]

         bootmom_temp = $
            measure_moments( x=bootx, y=booty, v=bootv, t=boott $
                       , /extrap, extarg=0, empty_props = empty_props, $
					   rms = rms, /do_bootstrap)
		 
		 this_mom2maj_extrap = bootmom_temp.mom2maj_extrap
		 this_mom2min_extrap = bootmom_temp.mom2min_extrap
		 this_momposang = bootmom_temp.momposang
		 this_mom2v_extrap = bootmom_temp.mom2v_extrap
		 this_mom2v_gauss_extrap = bootmom_temp.mom2v_gauss_extrap
		 this_mom0_extrap = bootmom_temp.mom0_extrap

		 ; DECONVOV BEAM
		 sig_to_fwhm = 2.354
 		 deconvolve_gauss $
           , meas_maj = this_mom2maj_extrap*sig_to_fwhm $  ; in pixel
           , meas_min = this_mom2min_extrap*sig_to_fwhm $  ; in pixel
           , meas_pa = this_momposang/!dtor $ ; INPUT AS DEGREES
           , beam_maj = bmaj_pix $
           , beam_min = bmin_pix $
           , beam_pa = bpa_deg $
           , src_maj = src_maj $
           , src_min = src_min $
           , src_pa = src_pa $
           , worked = worked $
           , point = point
        src_pa *= !dtor            ; SAVE AS RADIANS

		this_mom2maj_extrap_deconv = src_maj/sig_to_fwhm  ; in pixel
		this_mom2min_extrap_deconv = src_min/sig_to_fwhm  ; in pixel
		this_momposang_deconv = src_pa                    ; in radian

		; DECONVOV VELOCITY
		chantosig = 1.0/sqrt(2.0 * !pi)
		this_mom2v_extrap_deconv = sqrt(this_mom2v_extrap^2. - (1.0 * chantosig)^2.)
		this_mom2v_gauss_extrap_deconv = sqrt(this_mom2v_gauss_extrap^2. - (1.0 * chantosig)^2.)

		; RECORD NUMBERS
        mom2maj_extrap_deconv[j] = this_mom2maj_extrap_deconv
        mom2min_extrap_deconv[j] = this_mom2min_extrap_deconv
        momposang_deconv[j] = this_momposang_deconv
        mom2v_extrap_deconv[j] = this_mom2v_extrap_deconv
        mom2v_gauss_extrap_deconv[j] = this_mom2v_gauss_extrap_deconv 
        mom0_extrap[j] = this_mom0_extrap

	   endfor    ; end j = 0, ... bootstrap-1

	   mom2rad_extrap_deconv = sqrt(mom2maj_extrap_deconv * mom2min_extrap_deconv)
	   this_mom.radrms_extrap_deconv_uc = indfac * mad(mom2rad_extrap_deconv) $
		   /median(mom2rad_extrap_deconv)
	   this_mom.vrms_extrap_deconv_uc = indfac * mad(mom2v_extrap_deconv) $
		   /median(mom2v_extrap_deconv)
	   this_mom.vrms_gauss_extrap_deconv_uc = indfac * mad(mom2v_gauss_extrap_deconv) $
		   /median(mom2v_gauss_extrap_deconv)
	   this_mom.lum_extrap_uc = indfac * mad(mom0_extrap)/median(mom0_extrap)
	   this_mom.mass_extrap_uc = indfac * mad(mom0_extrap)/median(mom0_extrap)
       vircoeff = 1162.3 ; changed by Eric, was 1040.0
	   virmass_extrap_deconv = vircoeff * mom2rad_extrap_deconv * mom2v_extrap_deconv^2.
	   this_mom.virmass_extrap_deconv_uc = indfac * mad(virmass_extrap_deconv) $
		   /median(virmass_extrap_deconv)
	   virmass_gauss_extrap_deconv = vircoeff * mom2rad_extrap_deconv * mom2v_gauss_extrap_deconv^2.
	   this_mom.virmass_gauss_extrap_deconv_uc = indfac * mad(virmass_gauss_extrap_deconv) $
		   /median(virmass_gauss_extrap_deconv)

	 endif    ; end (bootstrap gt 0)  


     this_mom.peaknum = cloudlist[i]

     if n_elements(moments) eq 0 then begin
        moments = [this_mom]        
     endif else begin
        moments = [moments, this_mom]
     endelse

  endfor    ; end i = 0, ... nclouds-1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  save, file=outfile, moments

  ; mwrfits, moments, './test.fits', /create ; added by Eric

end                             ; OF CUBE_TO_MOMENTS
