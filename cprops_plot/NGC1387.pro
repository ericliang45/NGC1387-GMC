PRO NGC1387, steps = steps, asgn_mode=asgn_mode $
	           , noisebox = noisebox $
			   , hi_thresh = hi_thresh $
			   , hi_nchan = hi_nchan $
			   , lo_thresh = lo_thresh $
			   , lo_nchan = lo_nchan $
			   , min_pix = min_pix $
			   , min_area = min_area $
			   , minpix = minpix $
			   , minarea = minarea $
			   , minvchan = minvchan $
			   , snr = delta_is_snr $
			   , nonuniform = nonuniform $
			   , delta = delta $
			   , bootstrap = bootstrap $
			   , line_name = line_name $
			   , alpha = alpha $
			   , bclip = bclip $
			   , ingalpar = ingalpar $
			   , silent = silent

; Measure the GMCs properties of NGC1387
; EXAMPLE    NGC1387_GMC,steps=["MASK","LOCALMAX","ASSIGN","PROPS","OUTPUT"]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; 0. STEPS AND CLEAN
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ; 1) default steps and asgn_mode
  if n_elements(steps) eq 0 then steps=['OUT_GMC','OUT_LARSON1','OUT_ALPHA1'] ; ,'OUT_GMC','OUT_LARSON1','OUT_ALPHA1'
  if n_elements(asgn_mode) eq 0 then asgn_mode = ["CLFRIENDtoo"]
  ; 2) default noise cube parameters
  if n_elements(noisebox) eq 0 then noisebox = 3
  ; 3) default island/mask parameters
  if n_elements(hi_thresh) eq 0 then hi_thresh = 2
  if n_elements(hi_nchan) eq 0 then hi_nchan = 1
  if n_elements(lo_thresh) eq 0 then lo_thresh = 1.5
  if n_elements(lo_nchan) eq 0 then lo_nchan = 1
  if n_elements(min_pix) eq 0 then min_pix = 32    
  if n_elements(min_area) eq 0 then min_area = 32
  ; 4) default cloud decompositon parameters
  ;   pixel_per_beam = 16
  if n_elements(minpix) eq 0 then minpix = [160,144,128,112,96,80,64,48,32,16] ;[90,75,60,45,30,15]
  if n_elements(minarea) eq 0 then minarea = [160,144,128,112,96,80,64,48,32,16] ;[90,75,60,45,30,15] also have tested 1000, 500, 300, 200, 180, non-detection
  if n_elements(minvchan) eq 0 then minvchan = 2.0 
  if n_elements(minshape) eq 0 then minshape = 0.0 ; meaningless if smaller than minconvexity? was 0.45
  if n_elements(minconvexity) eq 0 then minconvexity = 0.50 ; tested, 0.40, 0.45, 0.50
  if n_elements(mintexture) eq 0 then mintexture = 0.0
  if n_elements(minUpar) eq 0 then minUpar = 0.00
  if n_elements(delta_is_snr) eq 0 then delta_is_snr = 0
  if n_elements(nonuniform) eq 0 then begin
	  if delta_is_snr eq 0 then nonuniform = 0
	  if delta_is_snr eq 1 then nonuniform = 1
  endif
  if n_elements(delta) eq 0 then delta = 2.0 ; in [K], 2*sigma (default), 1.5*sigma or 1*sigma; sigma of N1387 corrected cube = 1.95 K; uncorrected cube 1.1 K
  ;if n_elements(bclip) eq 0 then bclip = 1e4 ; correspond Tclip in [K]
  ; 5) default bootstrap parameters
  if n_elements(bootstrap) eq 0 then bootstrap = 1000 ;1000
  ; 6) default observation parameters
  if n_elements(line_name) eq 0 then line_name = "CO2-1"
  ; if n_elements(alpha) eq 0 then alpha = 4.8
  ; 7) default ingalpar file
  if n_elements(ingalpar) eq 0 then ingalpar = "./data/NGC1387_galpar.dat"

  ; bootstrap = 3, 10, 100, 1000; time (s) = 35, 50, 260, 2400; linear behavior, 2.3 sec per bootstrap

  ; "PREP", "NOISE", "MASK"
  ; "LOCALMAX","ASSIGN","PROPS" ("MOMENTS" & "PROPERTIES")
  ; "OUT_GMC", "OUT_PROPERTY_TABLE", "COMPARE_SIGMA", "OUT_MASS_FUNCTION", "OUT_ALPHA1", "OUT_LARSON1"
  ; "OUT_ANGMOM_COMPARE" "OUT_VFIELD" 

  number = ''

  ; for appendix B
  ; bootstrap = 100
  ; steps=["LOCALMAX","ASSIGN","PROPS","OUT_GMC"]

  ; minconvexity = 0.45
  ; delta = 3.0
  ; number = '/appendix/1'
  ; number_old = '_1'

  ; minconvexity = 0.45
  ; delta = 2.0
  ; number = '/appendix/2'
  ; number_old = '_2'

  ; minconvexity = 0.45
  ; delta = 1.0 ; was 4.0
  ; number = '/appendix/3'
  ; number_old = '_3'

  ; minconvexity = 0.55 ; was 0.40
  ; delta = 3.0
  ; number = '/appendix/4'
  ; number_old = '_4'

  ; minconvexity = 0.55 ; was 0.40
  ; delta = 2.0
  ; number = '/appendix/5'
  ; number_old = '_5'

  ; minconvexity = 0.55 ; was 0.40
  ; delta = 1.0 ; was 4.0
  ; number = '/appendix/6'
  ; number_old = '_6'

  ; minconvexity = 0.50
  ; delta = 3.0
  ; number = '/appendix/7'
  ; number_old = '_7'

  ; minconvexity = 0.50
  ; delta = 2.0
  ; number = '/appendix/8' ; the chosen one
  ; number_old = '_8'

  ; minconvexity = 0.50
  ; delta = 1.0 ; was 4.0
  ; number = '/appendix/9'
  ; number_old = '_9'

  FILE_MKDIR, 'measurements'+number
  FILE_MKDIR, 'output'+number

; Read galaxy parameters from ingalpar file
  pa=0.0
  inc = 0.0
  ra_c = ''
  dec_c = ''
  dist_pc = 0.0
  blank = ''
  vsys = 0.0
  inner=300.
  outer=600.
  alpha = 0.

  openr,lun,ingalpar,/get_lun
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,pa
  readf,lun,blank
  readf,lun,inc
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,ra_c
  readf,lun,blank
  readf,lun,dec_c
  readf,lun,blank
  readf,lun,dist_pc
  readf,lun,blank
  readf,lun,vsys
  readf,lun,blank
  readf,lun,alpha

  close,lun
  free_lun,lun

  print, 'alpha_co: ', alpha
  

; 0.0 set up steps
  check1 = where((steps eq "ALL") or (steps eq "CLEAN"), ct1)
  if ct1 gt 0 then begin
  ; 0.1 clean './*.eps'
    ;figures=file_search('./*.eps')
    ;if (size(figures))[0] gt 0 then begin
  	;  file_delete,figures
    ;endif	  
  
  ; 0.2 clean './measurements' and './output' folder
    if keyword_set(silent) then begin
		delete = 1
	endif else begin
       enter=''
       read, enter, prompt='Delete all files in folder ./measurements?'
       pastr = ['y','Y','Yes','YES','yes']
       delete = 0
       if total(strmatch(pastr, enter)) gt 0 then begin
            enter2=''
            read, enter2, prompt='Are you sure ?'
            if total(strmatch(pastr, enter2)) gt 0 then begin
                delete =  1
            endif
       endif
	endelse
  
    if (delete eq 1) then begin
      files=file_search('./measurements/*')
      if (size(files))[0] gt 0 then begin
        file_delete,files
      endif
  
      files=file_search('./output/*')
      if (size(files))[0] gt 0 then begin
        file_delete,files
      endif
    endif
  endif   ; END STEP "CLEAN"
  

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; II. PREPARATION
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; The cube needs to have km/s velocity units with units of Kelvin in
; order to proceed. "prep_cube" is a necessarily imperfect attempt to
; clean up a header to CPROP standards. 


; 2.1 CLEAN UP UNITS
  check2 = where((steps eq "ALL") or (steps eq "PREP"), ct2)
  if ct2 gt 0 then begin
    prep_cube $
       , in_file="./data/NGC1387_combine_clean5.image.pbcor-VEL.fits" $
       , out_file="./measurements"+number+"/NGC1387_CO21_cube_2kms_correct.fits" $
  	   , dist_pc= dist_pc $
	   , xco = alpha $ ; not stored anywhere, only for printing output
	   , galaxy_cen = [ra_c, dec_c] $
       , line_name= line_name
  
  ; 2.2 CHECK THE HEADER
    print, cprops_check_header(in_file="./measurements"+number+"/NGC1387_CO21_cube_2kms_correct.fits" $
                               , perfect=perfect $
                               , comments=comments)
  
  ; 2.3 PRINT ANY COMMENTS FROM THE HEADER CHECK
    print, comments
    print, 'pass'
  endif   ; END STEP "PREP"

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; III. GENERATE NOISE CUBE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; Estimate the noise from the data. We will use this to make a mask
; identifying bright emission. The noise estimate can have a wide
; range of complexities. Here the call uses a moving box of 9x9 pixels
; and calculates a two-d noise map. Optionally, it can fit a spectral
; dependence of the noise or just treat the whole cube as a single
; distribution. The example case here would probably be better served
; with a single value.

  check3 = where((steps eq "ALL") or (steps eq "NOISE"), ct3)
  if ct3 gt 0 then begin
    make_noise_cube $
       , cube_file = "./measurements"+number+"/NGC1387_CO21_cube_2kms_correct.fits" $
       , out_file = "./measurements"+number+"/NGC1387_CO21_cube_2kms_noise.fits" $
       , box = noisebox $
       , /iterate $
       , /twod_only $
       ;, /show $
       , /collapse 
  endif   ; END STEP "NOISE"

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; IV. IDENTIFY GALAXY REGION AND MAKE A MASK OF BRIGHT EMISSION
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; 4.1 Generate a mask for galaxy region only. Don't forget to check the figures 
;     - 'smooth_before_and_after.eps' and 'mom0_before_and_after_mask.eps' 
;     to make sure the maked region for galaxy is all right.

  check4 = where((steps eq "ALL") or (steps eq "MASK"), ct4)
  if ct4 gt 0 then begin
   ;make_galaxy_mask, infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits',$
   ;   rmsfile='./measurements/NGC1387_CO21_cube_2kms_noise.fits',$
   ;   outmask='./measurements/NGC1387_CO21_cube_2kms_galaxy_mask.fits',$
   ;   fbeam=[3.0], svel=25, $
   ;   hi_thresh=3.0, lo_thresh=1.5, chmin=2 
   ;make_galaxy_masktoo, $
   ;   infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits',$
   ;   max_dist= 1900., ingalpar = ingalpar, $
   ;   outdist='./measurements/NGC1387_CO21_cube_2kms_dist.fits', $
   ;   outmask='./measurements/NGC1387_CO21_cube_2kms_galaxy_mask.fits'
      
  ; 4.2 Make a mask of bright emission. Require 'hi_nchan' channels
  ; above signal-to-noise 'hi_thresh' and expand these regions 
  ; into a surface with at least 'lo_nchan' channels above signal-to
  ; -noise 'lo_thresh'
    make_cprops_mask $
       , infile = "./measurements"+number+"/NGC1387_CO21_cube_2kms_correct.fits" $
       , rmsfile='./measurements'+number+'/NGC1387_CO21_cube_2kms_noise.fits' $
       , inmask = "./data/NGC1387_combine_clean9.mask.fits" $

       , outfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_mask.fits" $
	     , outisland="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_island.fits" $
       , outmomzero="./measurements"+number+"/NGC1387_CO21_cube_2kms_mom0.fits" $
       , outmomone="./measurements"+number+"/NGC1387_CO21_cube_2kms_mom1.fits" $
       , hi_thresh = hi_thresh $
       , hi_nchan = hi_nchan $
       , lo_thresh = lo_thresh $
       , lo_nchan = lo_nchan $
       , min_pix = min_pix $               ; minimum pixels for island
       , min_area = min_area  $              ; minimum areas for island
	   , /nonuniform
      ; , /use_float ; added by Eric

   
  ; 4.2 Rerun the masking with the cube "inverted" (flipped to negative) in
  ; order to check the stringency of the mask. The number of pixels in
  ; the mask here give an estimate of the false positive rate.
  ;  make_cprops_mask $
  ;     , infile = "./measurements/NGC1387_CO21_cube_2kms_correct.fits" $
  ;     , outfile="./measurements/NGC1387_CO21_cube_2kms_falsepos_test.fits" $
  ;     , rmsfile='./measurements/NGC1387_CO21_cube_2kms_noise.fits' $
  ;     , hi_thresh = hi_thresh $
  ;     , hi_nchan = hi_nchan $
  ;     , lo_thresh = lo_thresh $
  ;     , lo_nchan = lo_nchan $
  ;     , min_pix = min_pix $               ; minimum pixels for island
  ;     , min_area = min_area  $             ; minimum areas for island
  ;     , /invert    $
  ;     , /nonuniform

  ; 4.3 Display Mask
  themask = readfits("./measurements"+number+"/NGC1387_CO21_cube_2kms_mask.fits", hdm)
  peak_map = max(readfits("./measurements"+number+"/NGC1387_CO21_cube_2kms_correct.fits", $
              hdr)*themask,dim=3,/nan)
  loadct, 33
  disp, peak_map, /sq


  endif  ; END STEP "MASK"

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; V. IDENTIFY LOCAL MAXIMA IN THE CUBE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; 5.1 Identify local maxima in the cube. The program slides a 
; sub-cube through the masked region and identifies significant local
; maxima. The results are saved in either a text file or IDL save file
; for further future analysis.
;
; This step allows a good deal of user input in what defines a
; significant local maximum. DELTA defines the contrast required above
; any contour shared with another maximum (this can be defined in
; intensity or SNR units). MINPIX and MINAREA define the area and
; volume required to be uniquely associated with an individual maximum
; for it to be valid. SPECFRIENDS and FRIENDS define the search box in
; pixel units. Think about this as having knobs to turn the contrast
; with the background, the search area, and the region required for a
; maximum to be valid.
; Manipulate the FIND_LOCAL_MAX call and iterate until happy.
  check5 = where((steps eq "ALL") or (steps eq "LOCALMAX"), ct5)
  if ct5 gt 0 then begin

    find_local_max_loop $
       , infile="./measurements_publication/NGC1387_CO21_cube_2kms_correct.fits" $
       , inmask ="./measurements_publication/NGC1387_CO21_cube_2kms_mask.fits" $
       ; , rmsfile='./measurements'+number+'/NGC1387_CO21_cube_2kms_noise.fits' $ ; if using delta as Kelvin, then no need for RMS file, also nonuniform = 0, delta_is_snr = 0

       , text_out="./measurements"+number+"/NGC1387_CO21_cube_2kms_lmax.txt" $
       , /verbose $
       , delta=delta $
	   ;, bclip=bclip $
       , snr=delta_is_snr $
       , minpix=minpix $                  
       , minarea=minarea $                
       , minshape = minshape $
       , minconvexity = minconvexity $
       , mintexture = mintexture $
       , minvchan=minvchan $
       , friends=1 $
       , specfriends=1 $
	   , /patchmode $
       ;, fscale=1 $
       ;, sigdiscont=1  $
  	   , nonuniform = nonuniform $
       , /no_deriv_decimate   

  ; 5.2 Visualize the selected peaks:
    ; themask = readfits("./measurements_publication/NGC1387_CO21_cube_2kms_mask.fits",hdm)
    ; peak_map = max(readfits("./measurements_publication/NGC1387_CO21_cube_2kms_correct.fits", $
    ;             hdr)*themask,dim=3,/nan)
    ; loadct, 33
    ; disp, peak_map, /sq
    ; readcol, "./measurements"+number+"/NGC1387_CO21_cube_2kms_lmax.txt", $
    ;          comment="#", x, y
    ; loadct, 0
    ; oplot, x, y, color=255, ps=1

  endif   ; END STEP "LOCALMAX"


; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; VI. ASSIGN VOXELS IN THE CUBE TO LOCAL MAXIMA
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Based on the list of local maxima and the mask, "color" the cube
; with assignments of voxels to discrete objects. The result is a cube
; with pixels corresponding to object number. The current propgram
; includes two methodologies to do this: the "CPROPS" way of assigning
; only emission uniquely associated with a maximum (defined as
; emission above any shared contour with another maximum) and the
; "CLUMPFIND" way of assigning emission within a shared contour to the
; nearest maximum on the basis of proximity.

  check6 = where((steps eq "ALL") or (steps eq "ASSIGN"), ct6)
  if ct6 gt 0 then begin

  ; 6.4 ROSOLOWSKY CLUMPFIND FRIEND-BY-FRIEND ASSIGNMENT
    asgn_ind = where((asgn_mode eq "CLFRIEND") or (asgn_mode eq "CLFRIENDtoo") or (asgn_mode eq "ALL"), asgn_ct) ; edited by Eric
    if (asgn_ct ge 1) then begin
      assign_clfriend $ 
       , infile="./measurements_publication/NGC1387_CO21_cube_2kms_correct.fits" $
       , inmask="./measurements_publication/NGC1387_CO21_cube_2kms_mask.fits" $
       , kernfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_lmax.txt" $

       , outfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriend.fits" $ ; once changed by Eric

  	   ; , rmsfile='./measurements'+number+'/NGC1387_CO21_cube_2kms_noise.fits' $
       , /diskmode $
       , gal_par = [inc, pa] $
       , nonuniform = nonuniform 
    endif

  ; 6.5 ROSOLOWSKY CLUMPFIND FRIEND-BY-FRIEND ASSIGNMENT
  ;     Make sure all resolved clouds have minimum convexity
    asgn_ind = where((asgn_mode eq "CLFRIENDtoo") or (asgn_mode eq "ALL") or (asgn_mode eq 'CHECK'), asgn_ct) ; once changed by Eric
    if (asgn_ct ge 1) then begin
      assign_check_convexity $ 
       , infile="./measurements_publication/NGC1387_CO21_cube_2kms_correct.fits" $
       , inmask="./measurements_publication/NGC1387_CO21_cube_2kms_mask.fits" $
       , rmsfile='./measurements_publication/NGC1387_CO21_cube_2kms_noise.fits' $
       , kernfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_lmax.txt" $
       , inassign="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriend.fits" $

       , outfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits" $
  	   , minarea = minarea[-1]  $
       , minpix = minpix[-1]  $
       , minshape = minshape $
       , minconvexity = minconvexity $
       , mintexture = mintexture $
       , minUpar = minUpar $
       , delta = delta $
       , snr = delta_is_snr $
       , nonuniform = nonuniform
    endif

  endif   ; END STEP "ASSIGN"

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; VII. MEASURE CLOUD MOMENTS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; Using the data cube (note that we now use the primary-beam corrected
; cube) and the assignment cubes generated in the last step to measure
; moments for each object in the assignment cube. Here moments is used
; loosely to refer to pixelwise properties (rather than properties
; with physical units, which we compute in the next step). 
;
; Additional analysis code could also be plugged in here to operate on
; the assignment+data combination.

  check7 = where((steps eq "ALL") or (steps eq "PROPS") or (steps eq "MOMENTS"), ct7)
  if ct7 gt 0 then begin

  ; 7.5 THE CLFRIENDtoo CASE
    asgn_ind = where((asgn_mode eq "CLFRIENDtoo") or (asgn_mode eq "ALL"), asgn_ct)
    if (asgn_ct ge 1) then begin
      cube_to_moments $
       , infile="./measurements_publication/NGC1387_CO21_cube_2kms_correct.fits" $
       , rmsfile='./measurements_publication/NGC1387_CO21_cube_2kms_noise.fits' $
       , inassign="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits" $
       
       , outfile="./measurements"+number+"/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl" $
  	   , bootstrap = bootstrap $
       , /verbose
    endif
  endif


  ; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ; VIII. CONVERT MOMENTS TO PROPERTIES
  ; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ; We now have a structure of moments for each object. The last step
  ; before science is to convert this into physical units. This step
  ; accepts moment measurements and uses a header and other assumptions
  ; (e.g., distance, conversion factor, etc.) to write a "properties"
  ; IDL structure. We run it once each for the two assignment cubes.
  
  check72 = where((steps eq "ALL") or (steps eq "PROPS") or (steps eq "PROPERTIES"), ct7)
  if ct7 gt 0 then begin

  ; 8.5 THE CLFRIENDtoo Case
    asgn_ind = where((asgn_mode eq "CLFRIENDtoo") or (asgn_mode eq "ALL"), asgn_ct)
    if (asgn_ct ge 1) then begin
      print, 'dist', dist_pc, 'alpha', alpha, 'pa', pa, 'inc', inc, 'ra_c',ra_c,'dec_c',dec_c
      moments_to_props $
       , hdrfile="./measurements/NGC1387_CO21_cube_2kms_correct.fits" $
       , infile="../work/measurements"+number_old+"/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl" $
       ; , inbootfile="./measurements_test"+number+"/NGC1387_CO21_cube_2kms_bootmomra_clfriendtoo.idl" $

       , idl_file="./measurements"+number+"/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl" $
       , fits_file="./measurements"+number+"/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits" $

       , dist= dist_pc $
  	   , alpha= alpha $     
      , pa = pa $
      , inc = inc $
      , ra_c = ra_c $
      , dec_c = dec_c $  
       , /verbose
    endif


  endif  ; END STEP "PROPS"
  

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; VIIII. GENERATE GMC TABLE AND FIGURE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; 9.1 GENERATE GMC TABLE
  ; GMC properties table

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_TABLE") or (steps eq "OUT_PROPERTY_TABLE"), ct8_1_1)
  if ct8_1_1 gt 0 then begin

    cpropstoo_gmc_table, $
              ; infile='./measurements'+number+'/NGC1387_CO21_cube_2kms_correct.fits', $
               ; inpmom='./measurements'+number+'/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl', $
               inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl', $
               kinematic_table = './output'+number+'/NGC1387_angmom_comparison-v11.csv',$
				ingalpar='./data/NGC1387_galpar.dat', $
          

               outtable='./output_publication/NGC1387_gmc_table.csv'
  endif  ;END STEP "OUT_TABLE"


  ; Pex table
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_TABLE") or (steps eq "OUT_Pext_TABLE"), ct8_1_2)
  if ct8_1_2 gt 0 then begin
    cpropstoo_cal_pext, intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , mge_file = './data/NGC1387_mge_gauss.csv'   $
            , incs = inc  $
            , dist_pc = dist_pc  $
			;, sigma_mass = 36.   $
			, infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
            , inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits' $
            , indist = "./measurements/NGC1387_CO21_cube_2kms_dist.fits" $
			, alpha = alpha $
            ;, inner_pc = inner_pc  $
            ;, ring_pc = [500., 900]  $
            ;, outer_pc = [900., 1e6]  $
			, out_vz2_distribution = './output/NGC1387_vz2_distribution.eps' $
			, out_ZR_ratio_distribution = './output/NGC1387_ZR_ratio_distribution.eps' $
			, out_zeta_distribution = './output/NGC1387_zeta_distribution.eps' $
			, out_stellar_gas_mass_comparison = './output/NGC1387_Mstar_vs_Mgas.eps' $
			, out_pin_pex = './output/NGC1387_Pin_Pex_comparison.eps' $
			, out_pext_sigmaGas = './output/NGC1387_sigmaGas_pext.eps' $
			, out_rhoStar_vturb = './output/NGC1387_vturb_rhoStar.eps' $
            , outtable ='./output/NGC1387_cloud_Pext.csv'
  endif

  ; Shear table
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_TABLE") or (steps eq "OUT_SHEAR_TABLE"), ct8_1_3)
  if ct8_1_3 gt 0 then begin
	cpropstoo_shear_table $
			, infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits'$
			;, inassign='./measurements/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
            , intable ='./output'+number+'/NGC1387_gmc_table.csv'  $
			, ingalpar='./data/NGC1387_galpar.dat' $
            , inpmom='./measurements'+number+'/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
            , rotdat ='./data/NGC1387_velocity_curve.dat'   $
            ;, inner_pc = 500  $
            ;, ring_pc = [500., 900]  $
            ;, outer_pc = [900., 1e6]  $
			, out_beta_model_fig = './output'+number+'/NGC1387_beta_distribution.eps' $
			, out_beta_obser_fig = './output'+number+'/NGC1387_beta_observed_distribution.eps' $
			, /logbeta $
			, out_shear_r = './output'+number+'/NGC1387_A_vs_r.eps'  $
			, outtable = './output'+number+'/NGC1387_shear_table.csv'  $
			, out_vrms_compare = './output'+number+'/NGC1387_vrms_vs_vpred.eps'
  endif  ;END STEP "OUT_SHEAR_TABLE"

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
    (steps eq "OUT_TABLE") or (steps eq "COMPARE_SIGMA"), ct8_1_4)
  if ct8_1_4 gt 0 then begin
  compare_sigma $
      , infile='./measurements'+number+'/NGC1387_CO21_cube_2kms_correct.fits'$
      ;, inassign='./measurements/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
            ; , intable ='./output'+number+'/NGC1387_gmc_table.csv'  $
            , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
            , ingalpar='./data/NGC1387_galpar.dat' $
            , rotdat ='./data/NGC1387_velocity_curve.dat'   $
            ;, inner_pc = 500  $
            ;, ring_pc = [500., 900]  $
            ;, outer_pc = [900., 1e6]  $
      ; , out_beta_model_fig = './output'+number+'/NGC1387_beta_distribution.eps' $
      ; , out_beta_obser_fig = './output'+number+'/NGC1387_beta_observed_distribution.eps' $
      ; , /logbeta $
      ; , out_shear_r = './output'+number+'/NGC1387_A_vs_r.eps'  $
      , outtable = './output'+number+'/NGC1387_sig_diff.csv'
      ; , out_vrms_compare = './output'+number+'/NGC1387_vrms_vs_vpred.eps'
  endif  ;END STEP "OUT_SHEAR_TABLE"

; 9.2 GENERATE GMC FIGURE
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_GMC"), ct8_2)
  if ct8_2 gt 0 then begin
    cpropstoo_gmc_figure, infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits', $
             inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits', $
             ; inpmom='./measurements'+number+'/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl', $
             inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl', $
             inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits', $
              checkfile_re = './measurements'+number+'/NGC1387-resolved_bool.txt', $
                checkfile_unre = './measurements'+number+'/NGC1387-unresolved_bool.txt', $

             outfile='./output'+number+'/NGC1387_gmc_figure.eps', $
             /plcloudNr, $
             ; particular=[1022] , $
			 ;shear_table='./output/NGC1387_shear_table.csv',  $
			 ;plellradius = [500,1000], $
			 ingalpar='./data/NGC1387_galpar.dat', $
             impad = [50, 50, 20, 80], $
		     scale_label = '200 pc'

  endif  ;END STEP "OUT_GMC"


  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
    (steps eq "OUT_GMC_APP"), ct8_22)
  if ct8_22 gt 0 then begin
    cpropstoo_gmc_figure_app, infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits', $
             inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits', $
             inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl', $
             inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits', $
              checkfile_re = './measurements'+number+'/NGC1387-resolved_bool.txt', $
                checkfile_unre = './measurements'+number+'/NGC1387-unresolved_bool.txt', $
             outfile='./output'+number+'/NGC1387_gmc_figure-app.eps', $
             ; /plcloudNr, $
             ; particular=[1022] , $
       ingalpar='./data/NGC1387_galpar.dat', $
             impad = [50, 50, 20, 80], $
         scale_label = '200 pc'
  endif  ;END STEP "OUT_GMC"

; 9.3 Generate a figure for "GMC mass function distribution"
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_MASS_FUNCTION"), ct8_3)
  if ct8_3 gt 0 then begin
    cpropstoo_mass_function   $
			, infile='./measurements'+number+'/NGC1387_CO21_cube_2kms_correct.fits'  $
			, inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl'  $
            ; , intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , minmass = 5.126695759569771 $
			, dist_pc = dist_pc $
			, xco = alpha $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
            , /notrunc $
                          , checkfile_re = './measurements'+number+'/NGC1387-resolved_bool.txt' $ 
                          , checkfile_unre = './measurements'+number+'/NGC1387-unresolved_bool.txt' $ 
            , outfile='./output'+number+'/NGC1387_mass_function.eps'
  endif  ; END STEP "OUT_MASS_FUNCTION"

; 9.4 Generate a figure for "Larson's Relation"
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_LARSON") or (steps eq "OUT_LARSON1"), ct8_4_1)
  if ct8_4_1 gt 0 then begin
    cpropstoo_larson_relation $ ; , intable='./output'+number+'/NGC1387_gmc_table.csv'
            , inprop = './measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
            ,  checkfile_re = './measurements/NGC1387-resolved_bool.txt' $
              ; , /nofitting $
            , outfile='./output'+number+'/NGC1387_larson_relation_obs.eps' $
              , alpha_co = alpha
  endif  ; END STEP "OUT_LARSON"


  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
    (steps eq "OUT_LARSON") or (steps eq "OUT_LARSON1_APP"), ct8_4_12)
  if ct8_4_12 gt 0 then begin
    cpropstoo_larson_relation_app $ ; , intable='./output'+number+'/NGC1387_gmc_table.csv'
            , inprop = './measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
            ,  checkfile_re = './measurements'+number+'/NGC1387-resolved_bool.txt' $
              ; , /nofitting $
            , outfile='./output'+number+'/NGC1387_larson_relation_obs.eps' $
              , alpha_co = alpha
  endif  ; END STEP "OUT_LARSON"


  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_LARSON") or (steps eq "OUT_LARSON2"), ct8_4_2)
  if ct8_4_2 gt 0 then begin
    cpropstoo_larson_relation, intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
			, vrmsmode = "turb" $
            , outfile='./output'+number+'/NGC1387_larson_relation_turb.eps'
  endif  ; END STEP "OUT_LARSON"

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_LARSON") or (steps eq "OUT_LARSON3"), ct8_4_3)
  if ct8_4_3 gt 0 then begin
    cpropstoo_larson_relation, intable='./output/NGC1387_gmc_table.csv' $
            ;, inner_pc = 500  $
            ;, ring_pc = [500., 900]  $
            ;, outer_pc = [900., 1e6]  $
			, vrmsmode = "eff" $
			, inc = inc $
			, shear_table='./output/NGC1387_shear_table.csv'  $
            , outfile='./output/NGC1387_larson_relation_eff.eps'
  endif  ; END STEP "OUT_LARSON"


; 9.5 GENERATE alpha FIGURE  
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_ALPHA") or (steps eq "OUT_ALPHA1"), ct8_5_1)
  if ct8_5_1 gt 0 then begin
    cpropstoo_alpha_figure $ 
    ; , intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
              , alpha_co = alpha $
            ,  checkfile_re = './measurements/NGC1387-resolved_bool.txt' $
            , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
			, plotval = [3.5, 8.5] $
        ; , /nofitting $
            ,outfile='./output/NGC1387_alpha_figure_obs.eps'
  endif  ; END STEP "OUT_ALPHA1"

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
    (steps eq "OUT_ALPHA") or (steps eq "OUT_ALPHA_APP"), ct8_5_12)
  if ct8_5_12 gt 0 then begin
    cpropstoo_alpha_figure_app $ 
    ; , intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
              , alpha_co = alpha $
            ,  checkfile_re = './measurements'+number+'/NGC1387-resolved_bool.txt' $
            , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
      , plotval = [3.5, 8.5] $
        ; , /nofitting $
            ,outfile='./output'+number+'/NGC1387_alpha_figure_obs-app.eps'
  endif  ; END STEP "OUT_ALPHA1"

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_ALPHA") or (steps eq "OUT_ALPHA2"), ct8_5_2)
  if ct8_5_2 gt 0 then begin
    cpropstoo_alpha_figure,intable='./output'+number+'/NGC1387_gmc_table.csv' $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
            , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
			, vrmsmode = 'turb' $
			, plotval = [3.5, 8.5] $
            , outfile='./output'+number+'/NGC1387_alpha_figure_turb.eps'
  endif  ; END STEP "OUT_ALPHA2"

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_ALPHA") or (steps eq "OUT_ALPHA3"), ct8_5_3)
  if ct8_5_3 gt 0 then begin
    cpropstoo_alpha_figure,intable='./output/NGC1387_gmc_table.csv' $
            ;, inner_pc = 500  $
            ;, ring_pc = [500., 900]  $
            ;, outer_pc = [900., 1e6]  $
            , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
			, vrmsmode = 'eff' $
			, plotval = [3.5, 8.5] $
			, shear_table='./output/NGC1387_shear_table.csv'  $
            , outfile='./output/NGC1387_alpha_figure_eff.eps'
  endif  ; END STEP "OUT_ALPHA2"


; 9.6 GMC KINEMATICS ANALYSIS
  ; 9.6.1 Generate a figure for "Molecular Cloud's Kinematics"
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_GMC_KINEMATIC") or (steps eq "OUT_CLOUD_KINEMATIC"), ct8_6_1)
  if ct8_6_1 gt 0 then begin
    cpropstoo_cloud_kinematics  $
                   , cloudnum=315 $
                   , fitmode='linear' $
  ;                 , cloudnum=25 $
  ;                 , fitmode='parabolic' $
                   , smfac=0.0 $ ; changed from 1.0 to 0.0 by Eric
                   , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
                   , inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
                   , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
                   , outfile='./output'+number+'/NGC1387_cloud_kinematics.eps'
  endif  ; END STEP "OUT_KINEMATICS"

  ; "test" for an update of kinematic properties for field map
  ; "test2-" for a new version with "/remove_boundary" and with previous update implemented
  ; 9.6.4 Compare observed GMC angular momentum with model
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
  (steps eq "OUT_GMC_KINEMATIC") or (steps eq "OUT_ANGMOM_COMPARE"), ct8_6_4)
  if ct8_6_4 gt 0 then begin
   cpropstoo_angmom_comparison $
                   , infile='./measurements'+number+'/NGC1387_CO21_cube_2kms_correct.fits' $
                   , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
                   , inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
           ; , intable = './output'+number+'/NGC1387_gmc_table.csv' $ ; obsolete, not used
           , rotfits='./data/NGC1387_rotation_vfield_2kms.fits' $ ;  './measurements/NGC1387_ROT_mom1_cube.fits' $
                      , checkfile_re = './measurements/NGC1387-resolved_bool.txt' $ 
                        , checkfile_unre = './measurements/NGC1387-unresolved_bool.txt' $
            ; , inner_pc = inner  $
            ; , ring_pc = [inner, outer]  $
            ; , outer_pc = [outer, 1e6]  $
           ; , /remove_boundary $ ; prefer not to use this mode
           , /fast $ ; essential
           ; , /plot_only $
           , outfile='./output'+number+'/NGC1387_angmom_comparison.eps' $
          , fits_file='./output'+number+'/NGC1387_angmom_comparison.fits'  $ 
           , outtable='./output'+number+'/NGC1387_angmom_comparison.csv' 
  endif  ; END STEP "OUT_ANGMOM_COMPARE"

  ; 9.6.2 Generate a figure for GMC vfield
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_GMC_KINEMATIC") or (steps eq "OUT_VFIELD"), ct8_6_2)
  if ct8_6_2 gt 0 then begin
    cpropstoo_gmc_vfield  $
                   ; , infile='./measurements'+number+'/NGC1387_CO21_cube_2kms_correct.fits' $
        				   ; , intable = './output/NGC1387_angmom_comparison.csv' $
                   , ingalpar='./data/NGC1387_galpar.dat' $
                   , rotfits='./data/NGC1387_rotation_vfield_2kms.fits' $
                   , inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits' $
                   , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
                   , intable = './output/NGC1387_angmom_comparison.csv' $
                      , checkfile_re = './measurements/NGC1387-resolved_bool.txt' $ 
                   ; , inpmom='./measurements'+number+'/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
                   ; , inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
                   , impar=[30, 4.0, 1/120., 12.] $
				   , impad = [28,42,0,85] $ ; [60, 75, 40, 120] $ ; used to be [55,55, 25, 90]; [ x1, x2, bottom, top]
				   , dist_range = [50, 950]    $
            , boundary = [inner,outer] $
				   , nvel=30 $
			       , /plcloudNr $
                   , outfile='./output'+number+'/NGC1387_gmc_vfield-cldnr.eps'
  endif  ; END STEP "OUT_VFIELD"

  ; 9.6.3 Generate output mom1 cube
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_GMC_KINEMATIC")or (steps eq "OUT_MOM_CUBE"), ct8_6_3)
  if ct8_6_3 gt 0 then begin
    cpropstoo_make_momcube  $
                   , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
                   , inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits' $
                   , rotfits='./data/NGC1387_rotation_vfield_2kms.fits' $
			       , ingalpar='./data/NGC1387_galpar.dat' $
                   , xco = alpha $
                   , outMom1_cube='./measurements/NGC1387_CO21_mom1_cube.fits' $
                   , outMom2_cube='./measurements/NGC1387_CO21_mom2_cube.fits' $
                   , outROvcube='./measurements/NGC1387_ROT_mom1_cube.fits' 
  endif  ; END STEP "OUT_MOM_CUBE"

 
  ; 9.6.5 Compare observed GMC angular momentum with model that gas 
  ;       moves following galaxy rotation
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	(steps eq "OUT_GMC_KINEMATIC") or (steps eq "OUT_ANGMOM_COMPAREII"), ct8_6_5)
  if ct8_6_5 gt 0 then begin
	 cpropstoo_angmom_comparisonToo $
                   , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
			       , ingalpar='./data/NGC1387_galpar.dat' $
                   , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
                   , inpmom='./measurements/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
				   , intable = './output/NGC1387_gmc_table.csv' $
				   , shear_table='./output/NGC1387_shear_table.csv'  $
                   ;, inner_pc = 500  $
                   ;, ring_pc = [500., 900]  $
                   ;, outer_pc = [900., 1e6]  $
				   , outfile='./output/NGC1387_angmom_comparisonII.eps'
  endif  ; END STEP "OUT_ANGMOM_COMPARE"

; 9.7 Generate a figure for show correlation between sigma_vR^-1/2 and surface 
;     density for extragalactic GMC populations.
  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_DYNAMICAL_STATE") or (steps eq "OUT_DYNAMICAL_STATE1"), ct8_7_1)
  if ct8_7_1 gt 0 then begin
    cpropstoo_dynamical_state  $
                   , intable='./output'+number+'/NGC1387_gmc_table.csv' $
                   , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl'  $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
                   , outfile='./output'+number+'/NGC1387_dynamical_state_obs.eps'
  endif  

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_DYNAMICAL_STATE") or (steps eq "OUT_DYNAMICAL_STATE2"), ct8_7_2)
  if ct8_7_2 gt 0 then begin
    cpropstoo_dynamical_state  $
                    , intable='./output'+number+'/NGC1387_gmc_table.csv' $
                    , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl'  $
            , inner_pc = inner  $
            , ring_pc = [inner, outer]  $
            , outer_pc = [outer, 1e6]  $
					, vrmsmode='turb' $
                    , outfile='./output'+number+'/NGC1387_dynamical_state_turb.eps'
  endif  

  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_DYNAMICAL_STATE") or (steps eq "OUT_DYNAMICAL_STATE3"), ct8_7_3)
  if ct8_7_3 gt 0 then begin
    cpropstoo_dynamical_state  $
                    , intable='./output/NGC1387_gmc_table.csv' $
                    , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl'  $
                    ;, inner_pc = 500  $
                    ;, ring_pc = [500., 900]  $
                    ;, outer_pc = [900., 1e6]  $
					, vrmsmode='eff' $
				    , shear_table='./output/NGC1387_shear_table.csv'  $
                    , outfile='./output/NGC1387_dynamical_state_eff.eps'
  endif  



; 9.9 Generate a figure for show the gmc geometry structure various with
;     radius (distance to galaxy center)
;  check8 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
;	  (steps eq "OUT_GMC_GEOMETRY_R"), ct8_9)
;  if ct8_9 gt 0 then begin
;    cpropstoo_gmc_geometry_vs_r  $
;                   , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
;                   , intable='./output/NGC1387_gmc_table.csv' $
;                   , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
;                   , inpmom='./measurements/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
;                   , mode=['dposang','ellipticity','gmcr','mlum','vrms'] $
;				   , r_range=[200, 300, 400]
;  endif

; 9.10 Generate a figure for show the mom0, mom1 and spectrum of 
;	   individual clouds
  check9 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_PLOTMOM"), ct8_10)
  if ct8_10 gt 0 then begin
    print,'here'
	cpropstoo_cloud_plotmom $
               , cloudnum = [100]  $
               , infile='./measurements/NGC1387_CO21_cube_2kms_mom1.fits' $
               , rotfile='./data/NGC1387_rotation_vfield_2kms.fits' $
			   ; , rotfile='./measurements/NGC1387_CO21_mom1_cube.fits' $
			   ; THE BELOW ROTFILE IS USED FOR /CALCSHEAR
               ;, rotfile='./data/NGC1387_velocity_curve.dat' $
			   ;, /calcShear $
               , inassign='./measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
               , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
               , inpmom='./measurements'+number+'/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
			   , /pl_ellipse $
			   ;, /plRs $
			   ;, shear_table='./output/NGC1387_shear_table.csv'  $
               , outfile='./output'+number+'/NGC1387_cloud_plotmom2.eps'
  endif

; 9.11 GENERATE GMC FIGURE
  check11 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_DEPROJECT_GMC"), ct8_11)
  if ct8_11 gt 0 then begin
    cpropstoo_deproject_gmc $
			   , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
               , intable='./output/NGC1387_gmc_table.csv' $
               ;, inner_pc = 500  $
               ;, ring_pc = [500., 900]  $
               ;, outer_pc = [900., 1e6]  $
			   , ingalpar='./data/NGC1387_galpar.dat' $
               , inassign='./measurements/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' $
               , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
			   , minor_degree = 125 $
               , inpmom='./measurements/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
			   , /no_deproject $
               , /plcloudNr $
			   , outfile='./output/NGC1387_deproject_gmc_figure.eps' $
			   , out_dpa='./output/NGC1387_dPa_distribution.eps' $
			   , out_ar='./output/NGC1387_axis_ratio_distribution.eps' $
               , impad = [10,10,2,2]

  endif  ;END STEP "OUT_DEPROJECT_GMC"

; 9.12 GENERATE VARIATION OF PROPERTIES AS DISTANCE
  check11 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
	  (steps eq "OUT_VARIATION_R"), ct8_12)
  if ct8_12 gt 0 then begin
    cpropstoo_variation_r $ ; , intable='./output'+number+'/NGC1387_gmc_table.csv'
                , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
                , inner_pc = inner  $
                , ring_pc = [inner, outer]  $
                , outer_pc = [outer, 1e6]  $
		        		, dist_binsize = 100 $
                , checkfile_re = './measurements/NGC1387-resolved_bool.txt' $
                , out_Rc_distribution = './output'+number+'/NGC1387-histogram-Rc.eps' $
                , out_Mlum_distribution = './output'+number+'/NGC1387-histogram-Mlum.eps' $
                , out_vrms_distribution = './output'+number+'/NGC1387-histogram-vrms.eps' $
                , out_sigmaGas_distribution = './output'+number+'/NGC1387-histogram-sigmaGas.eps' $
                , out_alphavir_distribution = './output'+number+'/NGC1387-histogram-alphavir.eps' $
                  , alpha_co = alpha

                ; , out_sigmaGas_vs_Rc = './output'+number+'/NGC1387_sigmaGas_vs_Rc.eps'  $
                ; , out_sigmaGas_vs_Mlum = './output'+number+'/NGC1387_sigmaGas_vs_Mlum.eps' 
  endif  ;END STEP "OUT_VARIATION_R"


; 9.13 CALCULATE ROCHE RADIUS
  check13 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
      (steps eq "OUT_ROCHE_RADIUS"), ct8_13)
  if ct8_13 gt 0 then begin
    cpropstoo_roche_radius, intable='./output/NGC1387_gmc_table.csv' $
          , rotdat ='./data/NGC1387_velocity_curve.dat'   $
          , rad_pc = [10, 1000]  $
		  , dist_pc = dist_pc $
          ;, mode = 'Binney2008' $
          ;, mode = 'Binney2008II' $
          ;, mode = 'Bertin2008' $
          , mode = 'Tan2000' $
          , out_roche_radius= './output/NGC1387_roche_radius.eps'

  endif  ;END STEP "OUT_ROCHE_RADIUS"

; 9.14 PLOT CLOUD SPECTRA
  check14 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
      (steps eq "OUT_CLOUD_SPECTRA"), ct8_14)
  if ct8_14 gt 0 then begin
    cpropstoo_cloud_spectra $
           , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
          ,  intable='./output'+number+'/NGC1387_gmc_table.csv' $
           , inassign="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits" $
           , cloudRange= [200,600] $
           ;, cloudnum= [164,179,125,135,206,124] $
           , outfile="./output"+number+"/NGC1387_cloud_spectra.eps"
  endif  ;END STEP "OUT_CLOUD_SPECTRA"

; 9.15  GENERATE MOMENT-0/1/2 FITS
  check15 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
      (steps eq "OUT_MOMENTS"), ct8_15)
  if ct8_15 gt 0 then begin
    cpropstoo_masked_moments $
           , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
           , inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits' $
           ; the below three input parameters are for cloud ellipse
           ;, /plcloud  $
           , inprop='./measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
           , inpmom='./measurements/NGC1387_CO21_cube_2kms_moments_clfriendtoo.idl' $
           , impad = [45,45, 15, 40] $
           , clipRMS = 1e6 $                ; in [km/s]
           , scale_label = '200 pc' $
           , vsys = vsys $
           ;, outmomzero = "./output/NGC1387_CO21_cube_2kms_mom0.fits" $
           ;, outmomone = "./output/NGC1387_CO21_cube_2kms_mom1.fits" $
           ;, outmomtwo = "./output/NGC1387_CO21_cube_2kms_mom2.fits" $
           , outzero = "./output/NGC1387_CO21_cube_2kms_mom0.eps"  $
           , outone = "./output/NGC1387_CO21_cube_2kms_mom1.eps" $
           , outtwo = "./output/NGC1387_CO21_cube_2kms_mom2.eps"
  endif  ;END STEP "OUT_CLOUD_SPECTRA"

; 9.17 Analysis stableness of galaxy disk and GMCs
  check17 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
      (steps eq "OUT_STABLE_ANALYSIS"), ct8_17)
  if ct8_17 gt 0 then begin
  cpropstoo_stable_analysis  $
            , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
			;, /onlydisc $
            , intable='./output/NGC1387_gmc_table.csv' $
            , inmask='./measurements/NGC1387_CO21_cube_2kms_mask.fits' $
            , inMom2='./measurements/NGC1387_CO21_mom2_cube.fits' $
            , incs = inc $
            , indist = "./measurements/NGC1387_CO21_cube_2kms_dist.fits" $
            , xco = alpha $
            , dist_pc = dist_pc $
            , dist_binsize = 30. $
			, maxdist = 1500 $
            , mge_file = './data/NGC1387_mge_gauss.csv'   $
            , rotdat ='./data/NGC1387_velocity_curve.dat'   $
			;, rotStddat = './data/NGC1387_omegaStd.dat' $
		    , out_timescale_r = './Output/NGC1387_timescale_vs_r.eps' $
            , out_Sigma_r = './Output/NGC1387_Sigma_gas_vs_r.eps' $
            , out_fgas_r = './Output/NGC1387_fgas_vs_r.eps' $
            , out_vrms_disk_r = './Output/NGC1387_vrms_disk_vs_r.eps' $
            , out_energy_r = './Output/NGC1387_energy_rate_vs_r.eps' $
            , out_rho_r = './Output/NGC1387_rho_gas_vs_r.eps' $
            , out_H_r = './Output/NGC1387_H_vs_r.eps' $
            , out_compare_scale = './Output/NGC1387_compare_scale.eps' $
            , out_Larson_test = './Output/NGC1387_Larson_test.eps' $
            , out_Q_vs_r = './Output/NGC1387_Q_vs_r.eps' $
            , out_Mc_vs_r = './Output/NGC1387_Mc_vs_r.eps' $
            ;, out_balance_dynamics = './Output/NGC1387_balance_dynamics.eps' $
            , out_compare_einject_evir = './Output/NGC1387_einject_vs_evir.eps' $
            , outtable = './Output/NGC1387_stable_analysis.csv'
  endif   ; END STEP "OUT_Stable_ANALYSIS"

; 9.18 PLOT CLOUD MOM0 MAP
  check18 = where((steps eq "ALL") or (steps eq "OUTPUT") or $
      (steps eq "OUT_CLOUD_MAP"), ct8_18)
  if ct8_18 gt 0 then begin
    cpropstoo_cloud_map $
           , infile='./measurements/NGC1387_CO21_cube_2kms_correct.fits' $
           , intable='./output'+number+'/NGC1387_gmc_table.csv' $
           , inassign="./measurements"+number+"/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits" $
           , inprop='./measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.idl' $
           , cloudRange= [400,600] $
           ; , cloudnum = [320] $
           , ncloud = 81 $
		   , minsize = 10 $
		   , /sm $
           , outfile="./output"+number+"/NGC1387_cloud_map.eps"
  endif  ;END STEP "OUT_CLOUD_MAP"

 
END
